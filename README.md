swift-statistics
================

Descriptive statistics in Swift, including an implementation of statistical accumulators.

Accumulators allow a sample of data to be gathered incrementally
in a compact form for a statistical calculation. 

Two protocols are defined: `Accumulator` and `PairAccumulator`.
The statistical accumulators conform to one of these protocols.
Sum, minimum, maximum, mean, sample variance, and population variance
are implemented, along with the Pearson correlation coefficient for data pairs.

The `value` property of an accumulator holds the statistical value
calculated by that accumulator or `nil` if an insufficient number
of sample points have been accumulated. Sample points are accumulated
using the `add(_: Double)` method for Accumulators or the
`add(x: Double, y: Double, weight: Double)` method for PairAccumulators.

An extension of the Sequence protocol for `Double` sequences
provides methods to calculate statistics on the sequence elements.
These protocol extensions act only on `Double` sequences,
but statistics for sequences of other numeric types can be
calculated either using a for..in loop or using this idiom:

    let variance = mySequence.lazy.map({Double($0)}).variance()
    let cc = myTupleSequence.lazy.map({ (Double($0.0), Double($0.1)) }).pearsonCorrelation()

The algorithms used for mean and variance follow the method of
Welford (1962). Relevant references (via Wikipedia):

Donald E. Knuth (1998). The Art of Computer Programming, Vol. 2:
Seminumerical Algorithms, 3rd edn., p. 232. Boston: Addison-Wesley.

B. P. Welford (1962)."Note on a method for calculating corrected sums
of squares and products". Technometrics 4(3):419–420.

Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1983). Algorithms
for Computing the Sample Variance: Analysis and Recommendations.
The American Statistician 37, 242-247.

Ling, Robert F. (1974). Comparison of Several Algorithms for Computing
Sample Means and Variances. Journal of the American Statistical
Association, Vol. 69, No. 348, 859-866.

Examples
--------

Using point-by-point accumulation:

    var minimum = Minimum()
    var maximum = Maximum()
    var mean = Mean()
    for i in 0..<100 {
        let x = Double(arc4_random_uniform(100))
        minimum.add(x)
        maximum.add(x)
        mean.add(x)
    }
    print("min: \(minimum.value!), max: \(maximum.value!), mean: \(mean.value!)")

Using the Sequence extensions:

    let sample = [0, 0.1, 1.0, 0.9]
    let ave = sample.mean()
    let stdev = sqrt(sample.sampleVariance()!)
    
Calculating a statistic for the elements of an [Int]:

    let sample = [0, 1, 10, 9]
    let average = sample.lazy.map({Double($0)}).mean()
    print(average!)

License
-------
MIT License

