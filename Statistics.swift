//
// Copyright © 2016 Charles Kissinger. All rights reserved.
// MIT License
//

// Statistical accumulators.
//
// Accumulators allow the data for statistical calculations to be
// gathered incrementally in a compact form.
//
// Sum, minimun, maximum, mean, sample variance, and population variance
// are implemented, along with the Pearson correlation coefficient for
// data pairs.
//
// The algorithms used for mean and variance follow the method
// of Welford (1962). Relevant references (via Wikipedia):
//
// Donald E. Knuth (1998). The Art of Computer Programming, Vol. 2:
// Seminumerical Algorithms, 3rd edn., p. 232. Boston: Addison-Wesley.
//
// B. P. Welford (1962)."Note on a method for calculating corrected sums
// of squares and products". Technometrics 4(3):419–420.
//
// Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1983). Algorithms
// for Computing the Sample Variance: Analysis and Recommendations.
// The American Statistician 37, 242-247.
//
// Ling, Robert F. (1974). Comparison of Several Algorithms for Computing
// Sample Means and Variances. Journal of the American Statistical
// Association, Vol. 69, No. 348, 859-866.

// The statistical protocol extensions for SequenceType act
// only on sequences of Doubles, but statistics for other
// numeric types can be calculated either using a for ... in loop
// or using this idiom:
//
//  let variance = mySequence.lazy.map({Double($0)}).variance()
//  let cc = mySequence.lazy.map({ (Double($0.0), Double($0.1)) }).pearsonCorrelation()

import Foundation // for sqrt()

/// An accumulator for calculating statistical values.
///
/// Accumulators allow incremental gathering of a sample of data
/// for the purpose of calculating a statistic describing the data.
public protocol AccumulatorBase: CustomStringConvertible,
    CustomDebugStringConvertible {

    /// The statistical value calculated by this accumulator
    /// or `nil` if an insufficient number of sample points
    /// have been accumulated.
    var value: Double? { get }

    /// The number of sample points accumulated.
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
        return "(value: \(value?.description ?? "nil") count: \(count))"
    }
}

public protocol Accumulator: AccumulatorBase {
    /// Add a sample point to the accumulator.
    ///
    /// - Parameter x: The sample point to be added.
    mutating func add(x: Double)
}

public protocol PairAccumulator: AccumulatorBase {
    /// Add a pair of sample points to the accumulator.
    ///
    /// - Parameter x: The first of a pair of sample points.
    /// - Parameter y: The second of a pair of sample points.
    /// - Parameter weight: The weight assigned to this pair of sample points
    ///   during calculation of the accumulator's value.
    mutating func add(x x: Double, y: Double, weight: Double)
}

/// An accumulator whose `value` property holds the sum
/// of the accumulated sample points.
public struct Sum: Accumulator {
    private var sum: Double = 0
    private var n: Double = 0
    /// The sum of the values stored by this accumulator
    /// or `nil` if no sample points have been accumulated.
    public var value: Double? { return n > 0 ? sum : nil }
    public var count: Int { return Int(n) }
    public mutating func add(x: Double) {
        sum += x
        n += 1
    }
    public mutating func reset() { sum = 0; n = 0 }
}

/// An accumulator whose `value` property holds the maximum
/// of the accumulated sample points.
public struct Maximum: Accumulator {
    private var max = -Double.infinity
    private var n = 0
    /// The maximum of the values stored by this accumulator
    /// or `nil` if no sample points have been accumulated.
    public var value: Double? { return n > 0 ? max : nil }
    public var count: Int { return Int(n) }
    public mutating func add(x: Double) {
        if x > max {
            max = x
        }
        n += 1
    }
    public mutating func reset() { max = 0; n = 0 }
}

/// An accumulator whose `value` property holds the minimum
/// of the accumulated sample points.
public struct Minimum: Accumulator {
    private var min = Double.infinity
    private var n = 0
    /// The minimum of the values stored by this accumulator
    /// or `nil` if no sample points have been accumulated.
    public var value: Double? { return n > 0 ? min : nil }
    public var count: Int { return Int(n) }
    public mutating func add(x: Double) {
        if x < min {
            min = x
        }
        n += 1
    }
    public mutating func reset() { min = 0; n = 0 }
}

/// An accumulator whose `value` property holds the arithmetic mean
/// of the accumulated sample points.
public struct Mean: Accumulator {
    private var mean: Double = 0
    private var n: Double = 0
    /// The arithmetic mean of the values stored by this accumulator
    /// or `nil` if no sample points have been accumulated.
    public var value: Double? { return n > 0 ? mean : nil  }
    public var count: Int { return Int(n) }
    public mutating func add(x: Double) {
        n += 1
        mean += (x - mean) / n
    }
    public mutating func reset() { mean = 0; n = 0 }
}

/// An accumulator whose `value` property holds the population variance
/// of the accumulated sample points.
public struct PopulationVariance: Accumulator {
    private var mean: Double = 0
    private var m2: Double = 0
    private var n: Double = 0
    /// The population variance of the values stored by this accumulator
    /// or `nil` if no sample points have been accumulated.
    public var value: Double? { return n > 0 ? m2 / n : nil }
    public var count: Int { return Int(n) }
    public mutating func add(x: Double) {
        n += 1
        let delta = x - mean
        mean += delta / n
        m2 += delta * (x - mean)
    }
    public mutating func reset() { mean = 0; m2 = 0; n = 0 }
}

/// An accumulator whose `value` property holds the sample variance
/// of the accumulated sample points.
public struct SampleVariance: Accumulator {
    private var mean: Double = 0
    private var m2: Double = 0
    private var n: Double = 0
    /// The sample variance of the values stored by this accumulator
    /// or `nil` if less than two data points have been accumulated.
    public var value: Double? { return n > 1 ? m2 / (n - 1) : nil }
    public var count: Int { return Int(n) }
    public mutating func add(x: Double) {
        n += 1
        let delta = x - mean
        mean += delta / n
        m2 += delta * (x - mean)
    }
    public mutating func reset() { mean = 0; m2 = 0; n = 0 }
}

/// An accumulator whose `value` property holds the Pearson sample
/// product-moment correlation coefficient (r) for the accumulated
/// pairs of sample points.
public struct PearsonCorrelation: PairAccumulator {
    private var meanX: Double = 0
    private var meanY: Double = 0
    private var n: Double = 0
    private var varX: Double = 0
    private var varY: Double = 0
    private var weightSum: Double = 0
    private var r: Double = 0
    /// The correlation coefficient of the set of data pairs stored by
    /// this accumulator or `nil` if less than two pairs of sample points
    /// have been accumulated.
    public var value: Double? { return n > 1 ? r / sqrt(varX * varY) : nil }
    public var count: Int { return Int(n) }

    public mutating func add(x x: Double, y: Double, weight: Double = 1) {
        n += 1
        let nextWeightSum = weightSum + weight
        let deltaX = x - meanX
        let scaleDeltaX = deltaX * weight / nextWeightSum
        varX += scaleDeltaX * deltaX * weightSum
        meanX += scaleDeltaX
        let deltaY = y - meanY
        let scaleDeltaY = deltaY * weight / nextWeightSum
        meanY += scaleDeltaY
        varY += scaleDeltaY * deltaY * weightSum
        weightSum = nextWeightSum
        r += deltaX * deltaY * (n - 1) * weight / n
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

public extension SequenceType where Generator.Element == Double {

    private func accumulate<T: Accumulator>(accumulator: T) -> Double? {
        var a = accumulator
        for x in self {
            a.add(x)
        }
        return a.value
    }

    /// Returns the sum of the elements in the sequence
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*length of sequence*)
    public func sum() -> Double? {
        return accumulate(Sum())
    }

    /// Returns the arithmetic mean of the elements in the sequence
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*length of sequence*)
    public func mean() -> Double? {
        return accumulate(Mean())
    }

    /// Returns the population variance of the elements in the sequence.
    /// or `nil` if the sequence is empty.
    ///
    /// - Complexity: O(*length of sequence*)
    public func populationVariance() -> Double? {
        return accumulate(PopulationVariance())
    }

    /// Returns the sample variance of the elements in the sequence.
    /// or `nil` if the sequence holds less than two elements.
    ///
    /// - Complexity: O(*length of sequence*)
    public func sampleVariance() -> Double? {
        return accumulate(SampleVariance())
    }
}

public extension SequenceType where Generator.Element == (Double, Double) {
    /// Returns the Pearson sample correlation coefficient
    /// of the tuple elements in the sequence.
    ///
    /// This method returns `nil` if the sequence holds less than two elements.
    /// - Complexity: O(*length of sequence*)
    public func pearsonCorrelation() -> Double? {
        var cc = PearsonCorrelation()
        for (x, y) in self {
            cc.add(x: x, y: y)
        }
        return cc.value
    }
}
