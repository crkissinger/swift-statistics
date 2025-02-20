//
// Copyright © 2016-2020 Charles Kissinger. All rights reserved.
// MIT License
//

/// Statistical functions and statistical accumulators.
///
/// Sum, minimun, maximum, arithmetic mean, geometric mean,
/// sample variance, and population variance are
/// implemented, along with the Pearson correlation coefficient
/// for data pairs.
///
/// Accumulators allow the data for statistical calculations to be
/// gathered incrementally in a compact form. In this way
/// statistical analysis is not limited to data stored in a single array.
/// Instead, data can be collected from arbitrary data structures,
/// streamed from disk or produced transiently,
//
// The algorithms used for mean and variance follow the method
// of Welford (1962).
// Relevant references (via Wikipedia):
//
// B. P. Welford (1962)."Note on a method for calculating corrected sums
// of squares and products". Technometrics 4(3):419–420.
//
// Donald E. Knuth (1998). The Art of Computer Programming, Vol. 2:
// Seminumerical Algorithms, 3rd edn., p. 232. Boston: Addison-Wesley.
//
// Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1983). Algorithms
// for Computing the Sample Variance: Analysis and Recommendations.
// The American Statistician 37, 242-247.
//
// Ling, Robert F. (1974). Comparison of Several Algorithms for Computing
// Sample Means and Variances. Journal of the American Statistical
// Association, Vol. 69, No. 348, 859-866.

import Foundation // for pow()

private typealias InternalFloatType = Double

/// An accumulator for calculating statistical values.
///
/// Accumulators allow incremental gathering of a sample of data
/// for the purpose of calculating a statistic describing the data.
public protocol AccumulatorBase: CustomStringConvertible,
    CustomDebugStringConvertible {

//    #if arch(x86_64)
//        associatedtype InternalFloatType = Float80
//    #else
//        associatedtype InternalFloatType = Double
//    #endif

    /// The statistical value calculated by this accumulator
    /// or `nil` if an insufficient number of data points
    /// have been accumulated.
    var value: Double? { get }

    /// The number of data points accumulated.
    var count: Int { get }

    /// Reinitialize the accumulator to an empty state.
    mutating func reset()

    var description: String { get }
    var debugDescription: String { get }
}

extension AccumulatorBase {
    public var description: String {
        return value?.description ?? "nil"
    }

    public var debugDescription: String {
        return "\(value?.debugDescription ?? "nil") (n=\(count))"
    }
}

public protocol Accumulator: AccumulatorBase {
    /// Add a data point to the accumulator.
    ///
    /// - Parameter x: The data point to be added.
    mutating func add<T>(_ x: T) where T: BinaryFloatingPoint
}

public protocol PairAccumulator: AccumulatorBase {
    /// Add a pair of data points to the accumulator.
    ///
    /// - Parameter x: The first of a pair of data points.
    /// - Parameter y: The second of a pair of data points.
    /// - Parameter weight: The weight assigned to this pair of data points
    ///   during calculation of the accumulator's value.
    mutating func add<T>(x: T, y: T, weight: T) where T: BinaryFloatingPoint
}

/// An accumulator whose `value` property holds the sum
/// of the accumulated data points.
public struct Sum: Accumulator {
    private var sum: InternalFloatType = 0
    private var n: Int = 0
    /// The sum of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated.
    public var value: Double? { return n > 0 ? Double(sum) : nil }
    public var count: Int { return n }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        sum += InternalFloatType(x)
        n += 1
    }
    public mutating func reset() { sum = 0; n = 0 }
}

/// An accumulator whose `value` property holds the maximum
/// of the accumulated data points.
public struct Maximum: Accumulator {
    private var max = -InternalFloatType.infinity
    private var n = 0
    /// The maximum of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated.
    public var value: Double? { return n > 0 ? Double(max) : nil }
    public var count: Int { return n }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        let value = InternalFloatType(x)
        if value > max {
            max = value
        }
        n += 1
    }
    public mutating func reset() { max = 0; n = 0 }
}

/// An accumulator whose `value` property holds the minimum
/// of the accumulated data points.
public struct Minimum: Accumulator {
    private var min = InternalFloatType.infinity
    private var n = 0
    /// The minimum of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated.
    public var value: Double? { return n > 0 ? Double(min) : nil }
    public var count: Int { return n }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        let value = InternalFloatType(x)
        if value < min {
            min = value
        }
        n += 1
    }
    public mutating func reset() { min = 0; n = 0 }
}

/// An accumulator whose `value` property holds the arithmetic mean
/// of the accumulated data points.
public struct Mean: Accumulator {
    private var mean: InternalFloatType = 0
    private var n: InternalFloatType = 0
    /// The arithmetic mean of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated.
    public var value: Double? { return n > 0 ? Double(mean) : nil  }
    public var count: Int { return Int(n) }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        n += 1
        mean += (InternalFloatType(x) - mean) / n
    }
    public mutating func reset() { mean = 0; n = 0 }
}

/// An accumulator whose `value` property holds the geometic mean
/// of the accumulated data points.
public struct GeometricMean: Accumulator {
    private var product: InternalFloatType = 1
    private var n: Int = 0
    private var invalidData = false
    /// The geometric mean of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated or the
    /// sample included.non-positive values.
    public var value: Double? {
        return n == 0 || invalidData ? nil : Double(pow(product, 1 / InternalFloatType(n)))
    }
    public var count: Int { return n }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        if x <= 0 {
            invalidData = true
        }
        product *= InternalFloatType(x)
        n += 1
    }
    public mutating func reset() { product = 1; n = 0; invalidData = false }
}

/// An accumulator whose `value` property holds the sample variance
/// of the accumulated data points.
public struct SampleStandardDeviation: Accumulator {
    private var mean: InternalFloatType = 0
    private var m2: InternalFloatType = 0
    private var n: InternalFloatType = 0
    /// The sample variance of the values stored by this accumulator
    /// or `nil` if less than two data points have been accumulated.
    public var value: Double? { return n > 1 ? Double((m2 / (n - 1)).squareRoot()) : nil }
    public var count: Int { return Int(n) }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        n += 1
        let delta = InternalFloatType(x) - mean
        mean += delta / n
        m2 += delta * (InternalFloatType(x) - mean)
    }
    public mutating func reset() { mean = 0; m2 = 0; n = 0 }
}

/// An accumulator whose `value` property holds the population variance
/// of the accumulated data points.
public struct PopulationStandardDeviation: Accumulator {
    private var mean: InternalFloatType = 0
    private var m2: InternalFloatType = 0
    private var n: InternalFloatType = 0
    /// The population variance of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated.
    public var value: Double? { return n > 0 ? Double((m2 / n).squareRoot()) : nil }
    public var count: Int { return Int(n) }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        n += 1
        let delta = InternalFloatType(x) - mean
        mean += delta / n
        m2 += delta * (InternalFloatType(x) - mean)
    }
    public mutating func reset() { mean = 0; m2 = 0; n = 0 }
}

/// An accumulator whose `value` property holds the sample variance
/// of the accumulated data points.
public struct SampleVariance: Accumulator {
    private var mean: InternalFloatType = 0
    private var m2: InternalFloatType = 0
    private var n: InternalFloatType = 0
    /// The sample variance of the values stored by this accumulator
    /// or `nil` if less than two data points have been accumulated.
    public var value: Double? { return n > 1 ? Double(m2 / (n - 1)) : nil }
    public var count: Int { return Int(n) }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        n += 1
        let delta = InternalFloatType(x) - mean
        mean += delta / n
        m2 += delta * (InternalFloatType(x) - mean)
    }
    public mutating func reset() { mean = 0; m2 = 0; n = 0 }
}

/// An accumulator whose `value` property holds the population variance
/// of the accumulated data points.
public struct PopulationVariance: Accumulator {
    private var mean: InternalFloatType = 0
    private var m2: InternalFloatType = 0
    private var n: InternalFloatType = 0
    /// The population variance of the values stored by this accumulator
    /// or `nil` if no data points have been accumulated.
    public var value: Double? { return n > 0 ? Double(m2 / n) : nil }
    public var count: Int { return Int(n) }
    public mutating func add<T>(_ x: T) where T: BinaryFloatingPoint {
        n += 1
        let delta = InternalFloatType(x) - mean
        mean += delta / n
        m2 += delta * (InternalFloatType(x) - mean)
    }
    public mutating func reset() { mean = 0; m2 = 0; n = 0 }
}
/// An accumulator whose `value` property holds the Pearson sample
/// product-moment correlation coefficient (r) for the accumulated
/// pairs of data points.
public struct PearsonCorrelation: PairAccumulator {
    private var meanX: InternalFloatType = 0
    private var meanY: InternalFloatType = 0
    private var n: InternalFloatType = 0
    private var varX: InternalFloatType = 0
    private var varY: InternalFloatType = 0
    private var weightSum: InternalFloatType = 0
    private var r: InternalFloatType = 0
    /// The correlation coefficient of the set of data pairs stored by
    /// this accumulator or `nil` if less than two pairs of data points
    /// have been accumulated.
    public var value: Double? { return n > 1 ? Double(r / (varX * varY).squareRoot()) : nil }
    public var count: Int { return Int(n) }

    public mutating func add<T>(x: T, y: T, weight: T = 1) where T: BinaryFloatingPoint {
        n += 1
        let nextWeightSum = weightSum + InternalFloatType(weight)
        let deltaX = InternalFloatType(x) - meanX
        let scaleDeltaX = deltaX * InternalFloatType(weight) / nextWeightSum
        varX += scaleDeltaX * deltaX * weightSum
        meanX += scaleDeltaX
        let deltaY = InternalFloatType(y) - meanY
        let scaleDeltaY = deltaY * InternalFloatType(weight) / nextWeightSum
        meanY += scaleDeltaY
        varY += scaleDeltaY * deltaY * weightSum
        weightSum = nextWeightSum
        r += deltaX * deltaY * (n - 1) * InternalFloatType(weight) / n
    }

    public mutating func reset() {
        meanX = 0
        meanY = 0
        varX = 0
        varY = 0
        weightSum = 0
        n = 0
    }
}

public extension Sequence where Element: BinaryFloatingPoint {

    /// Fills the accumulator with the contents of the `Self` and returns its value.
    /// - Returns: The `value` property of the accumulator, which can be `nil`.
    private func accumulate<T: Accumulator>(_ accumulator: T) -> Double? {
        return reduce(into: accumulator, { $0.add($1) }).value
    }

    /// Returns the sum of the elements in the sequence.
    ///
    /// - Returns: The sum of the elements in the sequence
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func sum() -> Double? {
        return accumulate(Sum())
    }

    /// Returns the arithmetic mean of the elements in the sequence.
    ///
    /// - Returns: The arithmetic mean of the elements in the sequence
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func mean() -> Double? {
        return accumulate(Mean())
    }

    /// Returns the geometic mean of the elements in the sequence.
    ///
    /// - Returns: The geometic mean of the elements in the sequence
    /// or `nil` if the sequence is empty or it contains non-positive elements.
    ///
    /// - Warning: The sequence must contain only positive numeric values.
    ///  If the sequence contains any non-positive values, this function returns `nil`.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func geometricMean() -> Double? {
        return accumulate(GeometricMean())
    }

    /// Returns the sample variance of the elements in the sequence.
    ///
    /// - Returns: The sample variance of the elements in the sequence.
    /// or `nil` if the sequence holds less than two elements.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func variance() -> Double? {
        return accumulate(SampleVariance())
    }

    /// Returns the population variance of the elements in the sequence.
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func populationVariance() -> Double? {
        return accumulate(PopulationVariance())
    }

    /// Returns the sample standard deviation of the elements in the sequence.
    ///
    /// - Returns: The sample standard deviation of the elements in the sequence.
    /// or `nil` if the sequence holds less than two elements.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func standardDeviation() -> Double? {
        return accumulate(SampleStandardDeviation())
    }

    /// Returns the population standard deviation of the elements in the sequence.
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func populationStandardDeviation() -> Double? {
        return accumulate(PopulationStandardDeviation())
    }
}

public extension Sequence where Element == (Double, Double) {
    /// Returns the Pearson sample correlation coefficient (*r*)  for the data pairs.
    ///
    /// - Returns: The Pearson sample correlation coefficient
    /// of the tuple elements in the sequence, or `nil` if the sequence contains
    /// less than two data pairs.
    ///
    /// - Complexity: O(*n*), where *n* is the length of the sequence.
    func pearsonCorrelation() -> Double? {
        var cc = PearsonCorrelation()
        for (x, y) in self {
            cc.add(x: x, y: y)
        }
        return cc.value
    }
}

/// Returns the Pearson sample correlation coefficient (*r*) for the corresponding data in two sequences.
///
/// If the two sequences are different lengths, the result is calculated over data pairs up to the end of the shorter sequence.
///
///  - Parameter sequence1: The first sequence of data points
///  - Parameter sequence2: The second sequence of data points
///
/// - Returns: The Pearson sample correlation coefficient
/// or `nil` if the sequences contains less than two data pairs.
///
/// - Complexity: O(*n*), where *n* is the length of the shorter sequence.
func pearsonCorrelation<T>( _ sequence1: T, _ sequence2: T) -> Double?
        where T: Sequence, T.Element: BinaryFloatingPoint {
    var cc = PearsonCorrelation()
    let zippedSequence = zip(sequence1, sequence2)
    for (x, y) in zippedSequence {
        cc.add(x: x, y: y)
    }
    return cc.value
}
