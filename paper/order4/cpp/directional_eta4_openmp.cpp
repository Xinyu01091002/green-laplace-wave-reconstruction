#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

using Complex = std::complex<double>;

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
Vec2 operator*(double value, Vec2 a) { return {value * a.x, value * a.y}; }
double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
double radial(Vec2 a) { return std::hypot(a.x, a.y); }

struct Mode2 {
    int32_t x = 0;
    int32_t y = 0;
};

Mode2 operator+(Mode2 a, Mode2 b) { return {a.x + b.x, a.y + b.y}; }
Mode2 operator-(Mode2 a, Mode2 b) { return {a.x - b.x, a.y - b.y}; }

struct Partition {
    std::array<int, 4> block{};
    int count = 0;
};

struct KernelOutput {
    Complex eta{};
    Complex psi{};
    Complex etaCoordinate{};
    Complex psiCoordinate{};
    Complex retained{};
    double normalizedDetuning = 0.0;
    double pairedResidual = 0.0;
    uint64_t strictZeroIntermediateCount = 0;
    double strictZeroSourceMaximum = 0.0;
    bool nearResonant = false;
};

struct DirectionalKernel {
    double gravity;
    double depth;
    double resonanceThreshold;
    std::vector<Partition> partitions[16][5];

    DirectionalKernel(double g, double h, double threshold)
        : gravity(g), depth(h), resonanceThreshold(threshold) {
        for (int mask = 1; mask < 16; ++mask) {
            std::vector<int> slots;
            for (int slot = 0; slot < 4; ++slot) {
                if ((mask & (1 << slot)) != 0) slots.push_back(slot);
            }
            for (int blockCount = 1; blockCount <= 4; ++blockCount) {
                int assignmentCount = 1;
                for (size_t index = 0; index < slots.size(); ++index) {
                    assignmentCount *= blockCount;
                }
                for (int code = 0; code < assignmentCount; ++code) {
                    int current = code;
                    int seen = 0;
                    Partition partition;
                    partition.count = blockCount;
                    for (int slot : slots) {
                        const int assignment = current % blockCount;
                        current /= blockCount;
                        seen |= 1 << assignment;
                        partition.block[assignment] |= 1 << slot;
                    }
                    if (seen == (1 << blockCount) - 1) {
                        partitions[mask][blockCount].push_back(partition);
                    }
                }
            }
        }
    }

    double q(Vec2 k) const {
        const double magnitude = radial(k);
        return magnitude * std::tanh(depth * magnitude);
    }

    double omega(Vec2 k) const { return std::sqrt(gravity * q(k)); }

    double zd(Vec2 k, int order) const {
        const double magnitude = radial(k);
        const double squared = magnitude * magnitude;
        switch (order) {
            case 0: return 1.0;
            case 1: return q(k);
            case 2: return squared;
            case 3: return squared * q(k);
            case 4: return squared * squared;
            default: throw std::runtime_error("unsupported vertical derivative order");
        }
    }

    std::array<Complex, 2> solve(
        const std::array<Complex, 2>& source, Vec2 k, double forcingOmega) const {
        const double qValue = q(k);
        const double denominator = forcingOmega * forcingOmega - gravity * qValue;
        return {
            (Complex(0.0, forcingOmega) * source[0] - qValue * source[1]) /
                denominator,
            (gravity * source[0] + Complex(0.0, forcingOmega) * source[1]) /
                denominator
        };
    }

    double g1(Vec2 kEta, Vec2 kPsi) const {
        return -q(kEta + kPsi) * q(kPsi) + dot(kPsi, kPsi) + dot(kEta, kPsi);
    }

    std::array<Complex, 2> quadraticNormalMap(
        Vec2 kA, int signA, Vec2 kB, int signB) const {
        const double qA = q(kA);
        const double qB = q(kB);
        const double omegaA = omega(kA);
        const double omegaB = omega(kB);
        const double mA = std::sqrt(omegaA / (2.0 * gravity));
        const double nA = std::sqrt(gravity / (2.0 * omegaA));
        const Complex etaA = mA;
        const Complex psiA(0.0, -signA * nA);
        const Complex etaB = 0.5;
        const Complex psiB(0.0, -signB * omegaB / (2.0 * qB));
        const std::array<Complex, 2> source{
            g1(kA, kB) * etaA * psiB + g1(kB, kA) * etaB * psiA,
            (dot(kA, kB) + qA * qB) * psiA * psiB
        };
        return solve(source, kA + kB, signA * omegaA + signB * omegaB);
    }

    std::array<Complex, 2> quadraticPhysicalCross(
        Vec2 kA, double omegaA, Complex etaA, Complex psiA,
        Vec2 kB, double omegaB, Complex etaB, Complex psiB) const {
        const std::array<Complex, 2> source{
            g1(kA, kB) * etaA * psiB + g1(kB, kA) * etaB * psiA,
            (dot(kA, kB) + q(kA) * q(kB)) * psiA * psiB
        };
        if (radial(kA + kB) == 0.0) return {};
        return solve(source, kA + kB, omegaA + omegaB);
    }

    KernelOutput evaluate(
        const std::array<Vec2, 4>& parentK,
        const std::array<int, 4>& signs) const {
        static constexpr double factorial[5] = {1.0, 1.0, 2.0, 6.0, 24.0};
        std::array<Vec2, 4> signedK{};
        std::array<double, 4> signedOmega{};
        for (int slot = 0; slot < 4; ++slot) {
            signedK[slot] = double(signs[slot]) * parentK[slot];
            signedOmega[slot] = signs[slot] * omega(parentK[slot]);
        }

        Vec2 subsetK[16]{};
        double subsetOmega[16]{};
        int degree[16]{};
        int charge[16]{};
        for (int mask = 1; mask < 16; ++mask) {
            for (int slot = 0; slot < 4; ++slot) {
                if ((mask & (1 << slot)) != 0) {
                    subsetK[mask] = subsetK[mask] + signedK[slot];
                    subsetOmega[mask] += signedOmega[slot];
                    degree[mask] += 1;
                    charge[mask] += signs[slot];
                }
            }
        }

        Complex eta[16]{};
        Complex phi[16]{};
        Complex psi[16]{};
        std::array<Complex, 2> surfaceHorizontal[16]{};
        Complex surfaceVertical[16]{};
        Complex retainedV3[16]{};
        std::array<Complex, 2> effectiveFourthSource{};
        std::array<Complex, 2> coordinateFourth{};
        uint64_t strictZeroIntermediateCount = 0;
        double strictZeroSourceMaximum = 0.0;

        const auto etaProduct = [&](const Partition& partition) {
            Complex value = 1.0;
            for (int block = 1; block < partition.count; ++block) {
                value *= eta[partition.block[block]];
            }
            return value;
        };

        const auto horizontalSurface = [&](int mask) {
            std::array<Complex, 2> value{};
            for (int order = 0; order < degree[mask]; ++order) {
                for (const auto& partition : partitions[mask][order + 1]) {
                    const int potentialMask = partition.block[0];
                    const double vertical = zd(subsetK[potentialMask], order);
                    const Complex common = phi[potentialMask] * etaProduct(partition) /
                        factorial[order];
                    value[0] += Complex(0.0, subsetK[potentialMask].x * vertical) * common;
                    value[1] += Complex(0.0, subsetK[potentialMask].y * vertical) * common;
                }
            }
            return value;
        };

        const auto verticalSurface = [&](int mask) {
            Complex value = 0.0;
            for (int order = 0; order < degree[mask]; ++order) {
                for (const auto& partition : partitions[mask][order + 1]) {
                    const int potentialMask = partition.block[0];
                    value += zd(subsetK[potentialMask], order + 1) *
                        phi[potentialMask] * etaProduct(partition) / factorial[order];
                }
            }
            return value;
        };

        const auto surfaceTaylor = [&](int mask) {
            Complex value = 0.0;
            for (int order = 1; order < degree[mask]; ++order) {
                for (const auto& partition : partitions[mask][order + 1]) {
                    const int potentialMask = partition.block[0];
                    value += zd(subsetK[potentialMask], order) * phi[potentialMask] *
                        etaProduct(partition) / factorial[order];
                }
            }
            return value;
        };

        const auto lowerVertical = [&](int mask) {
            Complex value = 0.0;
            for (int order = 1; order < degree[mask]; ++order) {
                for (const auto& partition : partitions[mask][order + 1]) {
                    const int potentialMask = partition.block[0];
                    value += zd(subsetK[potentialMask], order + 1) * phi[potentialMask] *
                        etaProduct(partition) / factorial[order];
                }
            }
            return value;
        };

        const auto lowerTime = [&](int mask) {
            Complex value = 0.0;
            for (int order = 1; order < degree[mask]; ++order) {
                for (const auto& partition : partitions[mask][order + 1]) {
                    const int potentialMask = partition.block[0];
                    value += Complex(0.0, -subsetOmega[potentialMask]) *
                        zd(subsetK[potentialMask], order) * phi[potentialMask] *
                        etaProduct(partition) / factorial[order];
                }
            }
            return value;
        };

        const auto ordinarySource = [&](int mask) {
            Complex horizontal = 0.0;
            Complex kinetic = 0.0;
            for (const auto& partition : partitions[mask][2]) {
                const int left = partition.block[0];
                const int right = partition.block[1];
                horizontal += (surfaceHorizontal[left][0] *
                    Complex(0.0, subsetK[right].x) + surfaceHorizontal[left][1] *
                    Complex(0.0, subsetK[right].y)) * eta[right];
                kinetic += 0.5 * (surfaceHorizontal[left][0] *
                    surfaceHorizontal[right][0] + surfaceHorizontal[left][1] *
                    surfaceHorizontal[right][1] + surfaceVertical[left] *
                    surfaceVertical[right]);
            }
            return std::array<Complex, 2>{
                horizontal - lowerVertical(mask), lowerTime(mask) + kinetic
            };
        };

        const auto dr2v3Push = [&](int mask) {
            std::array<Complex, 2> push{};
            for (int singleton = 0; singleton < 4; ++singleton) {
                if ((mask & (1 << singleton)) == 0) continue;
                const int cubicMask = mask & ~(1 << singleton);
                if (std::norm(retainedV3[cubicMask]) == 0.0) continue;
                const auto pair = quadraticNormalMap(
                    subsetK[cubicMask], charge[cubicMask] > 0 ? 1 : -1,
                    signedK[singleton], signs[singleton]);
                push[0] += retainedV3[cubicMask] * pair[0];
                push[1] += retainedV3[cubicMask] * pair[1];
            }
            return push;
        };

        const auto dr2y31Coordinate = [&](int mask) {
            std::array<Complex, 2> coordinate{};
            for (int singleton = 0; singleton < 4; ++singleton) {
                if ((mask & (1 << singleton)) == 0) continue;
                const int cubicMask = mask & ~(1 << singleton);
                const int cubicCharge = charge[cubicMask];
                if (std::abs(cubicCharge) != 1) continue;
                if (std::abs(eta[cubicMask]) + std::abs(psi[cubicMask]) == 0.0) continue;
                const double cubicFrequency = omega(subsetK[cubicMask]);
                const Complex alignedEta = 0.5 * (eta[cubicMask] +
                    Complex(0.0, cubicCharge * cubicFrequency / gravity) *
                    psi[cubicMask]);
                const Complex oppositeEta = 0.5 * (eta[cubicMask] -
                    Complex(0.0, cubicCharge * cubicFrequency / gravity) *
                    psi[cubicMask]);
                const Complex alignedPsi = Complex(0.0, -gravity * cubicCharge /
                    cubicFrequency) * alignedEta;
                const Complex oppositePsi = Complex(0.0, gravity * cubicCharge /
                    cubicFrequency) * oppositeEta;
                const double singletonQ = q(signedK[singleton]);
                const Complex singletonEta = 0.5;
                const Complex singletonPsi(0.0, -signs[singleton] *
                    omega(parentK[singleton]) / (2.0 * singletonQ));
                const auto aligned = quadraticPhysicalCross(
                    subsetK[cubicMask], cubicCharge * cubicFrequency,
                    alignedEta, alignedPsi, signedK[singleton], signedOmega[singleton],
                    singletonEta, singletonPsi);
                const auto opposite = quadraticPhysicalCross(
                    subsetK[cubicMask], -cubicCharge * cubicFrequency,
                    oppositeEta, oppositePsi, signedK[singleton], signedOmega[singleton],
                    singletonEta, singletonPsi);
                coordinate[0] += aligned[0] + opposite[0];
                coordinate[1] += aligned[1] + opposite[1];
            }
            return coordinate;
        };

        for (int currentDegree = 1; currentDegree <= 4; ++currentDegree) {
            for (int mask = 1; mask < 16; ++mask) {
                if (degree[mask] != currentDegree) continue;
                Complex taylor = 0.0;
                if (currentDegree == 1) {
                    eta[mask] = 0.5;
                    phi[mask] = Complex(0.0, -subsetOmega[mask] /
                        (2.0 * q(subsetK[mask])));
                    psi[mask] = phi[mask];
                } else {
                    const auto ordinary = ordinarySource(mask);
                    taylor = surfaceTaylor(mask);
                    std::array<Complex, 2> source{
                        -ordinary[0] - q(subsetK[mask]) * taylor,
                        -ordinary[1] - Complex(0.0, subsetOmega[mask]) * taylor
                    };
                    const double magnitude = radial(subsetK[mask]);
                    const bool strictZero = magnitude <= 100.0 *
                        std::numeric_limits<double>::epsilon() * std::max(magnitude, 1.0);
                    if (strictZero) {
                        eta[mask] = 0.0;
                        psi[mask] = 0.0;
                        strictZeroIntermediateCount += 1;
                        strictZeroSourceMaximum = std::max(strictZeroSourceMaximum,
                            std::hypot(std::abs(source[0]), std::abs(source[1])));
                    } else if (currentDegree == 2) {
                        const auto pair = solve(source, subsetK[mask], subsetOmega[mask]);
                        eta[mask] = pair[0];
                        psi[mask] = pair[1];
                    } else if (currentDegree == 3) {
                        if (std::abs(charge[mask]) == 3) {
                            const auto pair = solve(source, subsetK[mask], subsetOmega[mask]);
                            eta[mask] = pair[0];
                            psi[mask] = pair[1];
                        } else {
                            const double outputFrequency = omega(subsetK[mask]);
                            const double mScale = std::sqrt(outputFrequency /
                                (2.0 * gravity));
                            const double nScale = std::sqrt(gravity /
                                (2.0 * outputFrequency));
                            const int targetSign = charge[mask] > 0 ? 1 : -1;
                            retainedV3[mask] = source[0] / (2.0 * mScale) +
                                Complex(0.0, targetSign) * source[1] / (2.0 * nScale);
                            const Complex complement = source[0] / (2.0 * mScale) -
                                Complex(0.0, targetSign) * source[1] / (2.0 * nScale);
                            const double denominator = subsetOmega[mask] +
                                targetSign * outputFrequency;
                            eta[mask] = Complex(0.0, mScale) * complement / denominator;
                            psi[mask] = -targetSign * nScale * complement / denominator;
                        }
                    } else {
                        const auto push = dr2v3Push(mask);
                        effectiveFourthSource = {source[0] - push[0], source[1] - push[1]};
                        coordinateFourth = dr2y31Coordinate(mask);
                        const double qOutput = q(subsetK[mask]);
                        const double linear = gravity * qOutput;
                        const double temporal = subsetOmega[mask] * subsetOmega[mask];
                        const double rho = std::abs(linear - temporal) /
                            std::max(linear + temporal, std::numeric_limits<double>::min());
                        if (rho < resonanceThreshold) {
                            const double outputFrequency = omega(subsetK[mask]);
                            const double mScale = std::sqrt(outputFrequency /
                                (2.0 * gravity));
                            const double nScale = std::sqrt(gravity /
                                (2.0 * outputFrequency));
                            const int targetSign = subsetOmega[mask] >= 0.0 ? 1 : -1;
                            retainedV3[mask] = effectiveFourthSource[0] /
                                (2.0 * mScale) + Complex(0.0, targetSign) *
                                effectiveFourthSource[1] / (2.0 * nScale);
                            const Complex complement = effectiveFourthSource[0] /
                                (2.0 * mScale) - Complex(0.0, targetSign) *
                                effectiveFourthSource[1] / (2.0 * nScale);
                            const double denominator = subsetOmega[mask] +
                                targetSign * outputFrequency;
                            eta[mask] = Complex(0.0, mScale) * complement / denominator;
                            psi[mask] = -targetSign * nScale * complement / denominator;
                        } else {
                            const auto pair = solve(effectiveFourthSource,
                                subsetK[mask], subsetOmega[mask]);
                            eta[mask] = pair[0];
                            psi[mask] = pair[1];
                        }
                    }
                    phi[mask] = psi[mask] - taylor;
                }
                if (currentDegree < 4) {
                    surfaceHorizontal[mask] = horizontalSurface(mask);
                    surfaceVertical[mask] = verticalSurface(mask);
                }
            }
        }

        const double qOutput = q(subsetK[15]);
        const double linear = gravity * qOutput;
        const double temporal = subsetOmega[15] * subsetOmega[15];
        const double rho = std::abs(linear - temporal) /
            std::max(linear + temporal, std::numeric_limits<double>::min());
        Complex etaResidual = Complex(0.0, -subsetOmega[15]) * eta[15] -
            qOutput * psi[15] - effectiveFourthSource[0];
        Complex psiResidual = Complex(0.0, -subsetOmega[15]) * psi[15] +
            gravity * eta[15] - effectiveFourthSource[1];
        if (rho < resonanceThreshold) {
            const double outputFrequency = omega(subsetK[15]);
            const double mScale = std::sqrt(outputFrequency / (2.0 * gravity));
            const double nScale = std::sqrt(gravity / (2.0 * outputFrequency));
            const int targetSign = subsetOmega[15] >= 0.0 ? 1 : -1;
            etaResidual += mScale * retainedV3[15];
            psiResidual += Complex(0.0, -targetSign * nScale) * retainedV3[15];
        }
        return {
            eta[15], psi[15], coordinateFourth[0], coordinateFourth[1],
            retainedV3[15], rho,
            std::hypot(std::abs(etaResidual), std::abs(psiResidual)),
            strictZeroIntermediateCount, strictZeroSourceMaximum,
            rho < resonanceThreshold
        };
    }
};

template <typename T>
void readValue(std::ifstream& input, T& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!input) throw std::runtime_error("binary input read failed");
}

template <typename T>
void writeValue(std::ofstream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

struct ThreadAccumulation {
    std::vector<Complex> eta;
    std::vector<Complex> psi;
    std::vector<Complex> etaCoordinate;
    std::vector<Complex> psiCoordinate;
    std::vector<Complex> retained;
    explicit ThreadAccumulation(size_t size)
        : eta(size), psi(size), etaCoordinate(size), psiCoordinate(size),
          retained(size) {}
};

int repeatedFactor3(const std::array<uint16_t, 3>& indices) {
    if (indices[0] == indices[2]) return 6;
    if (indices[0] == indices[1] || indices[1] == indices[2]) return 2;
    return 1;
}

int repeatedFactor4(const std::array<uint16_t, 4>& indices) {
    int factor = 1;
    int run = 1;
    for (int index = 1; index < 4; ++index) {
        if (indices[index] == indices[index - 1]) {
            run += 1;
        } else {
            for (int value = 2; value <= run; ++value) factor *= value;
            run = 1;
        }
    }
    for (int value = 2; value <= run; ++value) factor *= value;
    return factor;
}

int main(int argc, char** argv) {
    try {
        if (argc < 3) {
            std::cerr << "input output [--sector 40|42|44]\n";
            return 2;
        }
        int sector = 40;
        for (int index = 3; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--sector" && index + 1 < argc) {
                sector = std::stoi(argv[++index]);
            } else if (argument.rfind("--sector=", 0) == 0) {
                sector = std::stoi(argument.substr(9));
            }
        }
        if (sector != 40 && sector != 42 && sector != 44) {
            throw std::runtime_error("sector must be 40, 42, or 44");
        }

        std::ifstream input(argv[1], std::ios::binary);
        char magic[8]{};
        input.read(magic, 8);
        if (std::string(magic, 8) != "DIR4IN01") {
            throw std::runtime_error("directional input magic mismatch");
        }
        uint64_t modeCount = 0;
        uint64_t nx = 0;
        uint64_t ny = 0;
        uint64_t caseCount = 0;
        readValue(input, modeCount);
        readValue(input, nx);
        readValue(input, ny);
        readValue(input, caseCount);
        double gravity = 0.0;
        double depth = 0.0;
        double deltaKx = 0.0;
        double deltaKy = 0.0;
        double threshold = 0.0;
        readValue(input, gravity);
        readValue(input, depth);
        readValue(input, deltaKx);
        readValue(input, deltaKy);
        readValue(input, threshold);
        std::vector<Mode2> modes(modeCount);
        for (auto& mode : modes) {
            readValue(input, mode.x);
            readValue(input, mode.y);
        }
        std::vector<std::vector<Complex>> amplitudes(
            caseCount, std::vector<Complex>(modeCount));
        for (auto& oneCase : amplitudes) {
            for (auto& value : oneCase) {
                double real = 0.0;
                double imaginary = 0.0;
                readValue(input, real);
                readValue(input, imaginary);
                value = {real, imaginary};
            }
        }

        const auto physicalK = [&](uint16_t index) {
            return Vec2{modes[index].x * deltaKx, modes[index].y * deltaKy};
        };
        const auto gridInside = [&](Mode2 mode) {
            return std::abs(mode.x) < static_cast<int32_t>(nx / 2) &&
                std::abs(mode.y) < static_cast<int32_t>(ny / 2);
        };
        const auto analyticHalfPlane = [](Mode2 mode) {
            return mode.x > 0 || (mode.x == 0 && mode.y > 0);
        };
        const auto location = [&](Mode2 mode, uint64_t caseIndex) {
            const uint64_t ix = static_cast<uint64_t>((mode.x % static_cast<int32_t>(nx) +
                static_cast<int32_t>(nx)) % static_cast<int32_t>(nx));
            const uint64_t iy = static_cast<uint64_t>((mode.y % static_cast<int32_t>(ny) +
                static_cast<int32_t>(ny)) % static_cast<int32_t>(ny));
            return caseIndex * nx * ny + ix + nx * iy;
        };

        const size_t stride = static_cast<size_t>(nx * ny * caseCount);
        const int threadCount = omp_get_max_threads();
        std::vector<ThreadAccumulation> threadOutput;
        threadOutput.reserve(threadCount);
        for (int thread = 0; thread < threadCount; ++thread) {
            threadOutput.emplace_back(stride);
        }
        DirectionalKernel kernel(gravity, depth, threshold);
        uint64_t totalCandidates = 0;
        uint64_t kept = 0;
        uint64_t near = 0;
        uint64_t strictZeroOutput = 0;
        uint64_t strictZeroIntermediate = 0;
        double maximumResidual = 0.0;
        double minimumDetuning = std::numeric_limits<double>::infinity();
        double maximumStrictZeroSource = 0.0;
        const double start = omp_get_wtime();

        const auto accumulate = [&](ThreadAccumulation& output, size_t out,
                                    Complex weight, const KernelOutput& value) {
            output.eta[out] += weight * value.eta;
            output.psi[out] += weight * value.psi;
            output.etaCoordinate[out] += weight * value.etaCoordinate;
            output.psiCoordinate[out] += weight * value.psiCoordinate;
            if (value.nearResonant) output.retained[out] += weight * value.retained;
        };

        if (sector == 40) {
            std::vector<std::array<uint16_t, 2>> pairs;
            for (uint16_t first = 0; first < modeCount; ++first) {
                for (uint16_t second = first; second < modeCount; ++second) {
                    pairs.push_back({first, second});
                }
            }
            totalCandidates = pairs.size() * pairs.size();
#pragma omp parallel for schedule(static, 1) reduction(+:kept,near,strictZeroOutput,strictZeroIntermediate) reduction(max:maximumResidual,maximumStrictZeroSource) reduction(min:minimumDetuning)
            for (int64_t minusIndex = 0;
                 minusIndex < static_cast<int64_t>(pairs.size()); ++minusIndex) {
                const int thread = omp_get_thread_num();
                const auto minus = pairs[minusIndex];
                for (const auto plus : pairs) {
                    const Mode2 target = modes[plus[0]] + modes[plus[1]] -
                        modes[minus[0]] - modes[minus[1]];
                    if (target.x == 0 && target.y == 0) {
                        strictZeroOutput += 1;
                        continue;
                    }
                    if (!analyticHalfPlane(target) || !gridInside(target)) continue;
                    const auto value = kernel.evaluate(
                        {physicalK(plus[0]), physicalK(plus[1]),
                         physicalK(minus[0]), physicalK(minus[1])},
                        {1, 1, -1, -1});
                    const int repeat = (1 + (plus[0] == plus[1])) *
                        (1 + (minus[0] == minus[1]));
                    for (uint64_t caseIndex = 0; caseIndex < caseCount; ++caseIndex) {
                        const auto& a = amplitudes[caseIndex];
                        const Complex weight = 2.0 * a[plus[0]] * a[plus[1]] *
                            std::conj(a[minus[0]]) * std::conj(a[minus[1]]) /
                            double(repeat);
                        accumulate(threadOutput[thread], location(target, caseIndex),
                            weight, value);
                    }
                    kept += 1;
                    near += value.nearResonant ? 1 : 0;
                    strictZeroIntermediate += value.strictZeroIntermediateCount;
                    maximumResidual = std::max(maximumResidual, value.pairedResidual);
                    minimumDetuning = std::min(minimumDetuning, value.normalizedDetuning);
                    maximumStrictZeroSource = std::max(maximumStrictZeroSource,
                        value.strictZeroSourceMaximum);
                }
            }
        } else if (sector == 42) {
            std::vector<std::array<uint16_t, 3>> triples;
            for (uint16_t first = 0; first < modeCount; ++first) {
                for (uint16_t second = first; second < modeCount; ++second) {
                    for (uint16_t third = second; third < modeCount; ++third) {
                        triples.push_back({first, second, third});
                    }
                }
            }
            totalCandidates = triples.size() * modeCount;
#pragma omp parallel for schedule(static, 1) reduction(+:kept,near,strictZeroOutput,strictZeroIntermediate) reduction(max:maximumResidual,maximumStrictZeroSource) reduction(min:minimumDetuning)
            for (int64_t minus = 0; minus < static_cast<int64_t>(modeCount); ++minus) {
                const int thread = omp_get_thread_num();
                for (const auto plus : triples) {
                    const Mode2 target = modes[plus[0]] + modes[plus[1]] +
                        modes[plus[2]] - modes[minus];
                    if (target.x == 0 && target.y == 0) {
                        strictZeroOutput += 1;
                        continue;
                    }
                    if (!gridInside(target)) continue;
                    const auto value = kernel.evaluate(
                        {physicalK(plus[0]), physicalK(plus[1]),
                         physicalK(plus[2]), physicalK(static_cast<uint16_t>(minus))},
                        {1, 1, 1, -1});
                    const int repeat = repeatedFactor3(plus);
                    for (uint64_t caseIndex = 0; caseIndex < caseCount; ++caseIndex) {
                        const auto& a = amplitudes[caseIndex];
                        const Complex weight = 2.0 * a[plus[0]] * a[plus[1]] *
                            a[plus[2]] * std::conj(a[minus]) / double(repeat);
                        accumulate(threadOutput[thread], location(target, caseIndex),
                            weight, value);
                    }
                    kept += 1;
                    near += value.nearResonant ? 1 : 0;
                    strictZeroIntermediate += value.strictZeroIntermediateCount;
                    maximumResidual = std::max(maximumResidual, value.pairedResidual);
                    minimumDetuning = std::min(minimumDetuning, value.normalizedDetuning);
                    maximumStrictZeroSource = std::max(maximumStrictZeroSource,
                        value.strictZeroSourceMaximum);
                }
            }
        } else {
            std::vector<std::array<uint16_t, 4>> quadruples;
            for (uint16_t first = 0; first < modeCount; ++first) {
                for (uint16_t second = first; second < modeCount; ++second) {
                    for (uint16_t third = second; third < modeCount; ++third) {
                        for (uint16_t fourth = third; fourth < modeCount; ++fourth) {
                            quadruples.push_back({first, second, third, fourth});
                        }
                    }
                }
            }
            totalCandidates = quadruples.size();
#pragma omp parallel for schedule(static) reduction(+:kept,near,strictZeroOutput,strictZeroIntermediate) reduction(max:maximumResidual,maximumStrictZeroSource) reduction(min:minimumDetuning)
            for (int64_t index = 0;
                 index < static_cast<int64_t>(quadruples.size()); ++index) {
                const int thread = omp_get_thread_num();
                const auto plus = quadruples[index];
                const Mode2 target = modes[plus[0]] + modes[plus[1]] +
                    modes[plus[2]] + modes[plus[3]];
                if (target.x == 0 && target.y == 0) {
                    strictZeroOutput += 1;
                    continue;
                }
                if (!gridInside(target)) continue;
                const auto value = kernel.evaluate(
                    {physicalK(plus[0]), physicalK(plus[1]),
                     physicalK(plus[2]), physicalK(plus[3])},
                    {1, 1, 1, 1});
                const int repeat = repeatedFactor4(plus);
                for (uint64_t caseIndex = 0; caseIndex < caseCount; ++caseIndex) {
                    const auto& a = amplitudes[caseIndex];
                    const Complex weight = 2.0 * a[plus[0]] * a[plus[1]] *
                        a[plus[2]] * a[plus[3]] / double(repeat);
                    accumulate(threadOutput[thread], location(target, caseIndex),
                        weight, value);
                }
                kept += 1;
                near += value.nearResonant ? 1 : 0;
                strictZeroIntermediate += value.strictZeroIntermediateCount;
                maximumResidual = std::max(maximumResidual, value.pairedResidual);
                minimumDetuning = std::min(minimumDetuning, value.normalizedDetuning);
                maximumStrictZeroSource = std::max(maximumStrictZeroSource,
                    value.strictZeroSourceMaximum);
            }
        }

        ThreadAccumulation merged(stride);
        for (const auto& oneThread : threadOutput) {
            for (size_t index = 0; index < stride; ++index) {
                merged.eta[index] += oneThread.eta[index];
                merged.psi[index] += oneThread.psi[index];
                merged.etaCoordinate[index] += oneThread.etaCoordinate[index];
                merged.psiCoordinate[index] += oneThread.psiCoordinate[index];
                merged.retained[index] += oneThread.retained[index];
            }
        }

        std::ofstream output(argv[2], std::ios::binary);
        const std::string outputMagic = sector == 40 ? "DIR40O01" :
            (sector == 42 ? "DIR42O01" : "DIR44O01");
        output.write(outputMagic.data(), 8);
        writeValue(output, nx);
        writeValue(output, ny);
        writeValue(output, caseCount);
        writeValue(output, totalCandidates);
        writeValue(output, kept);
        writeValue(output, near);
        writeValue(output, strictZeroOutput);
        writeValue(output, strictZeroIntermediate);
        const double seconds = omp_get_wtime() - start;
        writeValue(output, seconds);
        writeValue(output, maximumResidual);
        writeValue(output, minimumDetuning);
        writeValue(output, maximumStrictZeroSource);
        for (const auto* array : {&merged.eta, &merged.psi,
                                  &merged.etaCoordinate, &merged.psiCoordinate,
                                  &merged.retained}) {
            for (const Complex value : *array) {
                const double real = value.real();
                const double imaginary = value.imag();
                writeValue(output, real);
                writeValue(output, imaginary);
            }
        }
        std::cerr << "DIRECTIONAL_CPP_DONE sector=" << sector
                  << " cases=" << caseCount
                  << " total=" << totalCandidates
                  << " kept=" << kept
                  << " near=" << near
                  << " seconds=" << seconds
                  << " residual=" << maximumResidual
                  << " threads=" << threadCount << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DIRECTIONAL_CPP_ERROR " << error.what() << "\n";
        return 1;
    }
}
