import org.apache.commons.math3.distribution.UniformRealDistribution
import org.apache.commons.math3.distribution.GammaDistribution
import org.apache.commons.math3.distribution.ExponentialDistribution
import kotlin.math.roundToInt
// import kotlin.random.Random

const val synthesisSup35: Int = 700
const val synthesisHsp104: Int = 50
const val totalSynthesisRate = synthesisHsp104 + synthesisSup35

const val rho: Double = 31.03
const val shapeDaughter: Double = rho * 0.21 * 60
const val shape: Double = rho * 1.16 * 60
const val sizeThreshold: Int = 30
const val minAggSize: Int = 5

fun divTimeMother(): Double {
    return GammaDistribution(shape, 1 / rho).sample()
}

fun divTimeDaughter(): Double {
    return GammaDistribution(shapeDaughter, 1 / rho).sample() + divTimeMother()
}

class YeastCell(var sup35MoleculeCount: Int = 46680,
                var hsp104MoleculeCount: Int = 3326,
                var divisionCount: Int = 0,
                var nextDivisionTime: Double = divTimeDaughter(),
                var elapsedTime: Double = 0.0,
                var aggregateSizes: MutableList<Int> = mutableListOf(20, 20, 20, 40, 40, 40, 40, 50)
) {
    var hasPSI = (aggregateSizes.size > 0)
}

fun YeastCell.updateProteinLevel(timePoint:  Double) {
    val timePeriod = timePoint - this.elapsedTime
    this.sup35MoleculeCount += (timePeriod * synthesisSup35).roundToInt()
    this.hsp104MoleculeCount += (timePeriod * synthesisHsp104).roundToInt()
    this.elapsedTime = timePoint
}

fun YeastCell.divide(): YeastCell {

    // Creating a daughter cell object, with empty aggregate vector
    this.divisionCount++
    val daughterCell = YeastCell(aggregateSizes = mutableListOf<Int>(), elapsedTime = this.elapsedTime + 1)
    daughterCell.nextDivisionTime += this.elapsedTime

    // Splitting the aggregates between cells in accordance with size and 60/40 rule
    val remainingAggregates = mutableListOf<Int>()
    for (i in this.aggregateSizes) {
        val testValue = UniformRealDistribution(0.0, 1.0).sample()
        if (testValue < 0.4 && i < sizeThreshold) {
            daughterCell.aggregateSizes.add(i)
        } else {
            remainingAggregates.add(i)
        }
    }
    if (daughterCell.aggregateSizes.size > 0) {
        true.also { daughterCell.hasPSI = it }
    }
    if (remainingAggregates.size == 0) {
        false.also { this.hasPSI = it }
    }
    this.aggregateSizes = remainingAggregates

    daughterCell.sup35MoleculeCount = (0.4 * this.sup35MoleculeCount).roundToInt()
    daughterCell.hsp104MoleculeCount = (0.4 * this.hsp104MoleculeCount).roundToInt()

    daughterCell.sup35MoleculeCount = (0.6 * this.sup35MoleculeCount).roundToInt()
    daughterCell.hsp104MoleculeCount = (0.6 * this.hsp104MoleculeCount).roundToInt()

    this.nextDivisionTime += divTimeMother()
    this.elapsedTime += 1

    return daughterCell
}

fun YeastCell.nextReaction(conversionBeta: Double, fragmentationGamma: Double) {

    // Calculating reaction rates
    val conversionRate = conversionBeta * this.sup35MoleculeCount * this.aggregateSizes.size
    val fragmentSites = this.aggregateSizes.sum() - this.aggregateSizes.size
    val fragmentationRate = fragmentSites * fragmentationGamma * this.hsp104MoleculeCount / (this.hsp104MoleculeCount / 2 + fragmentSites)
    val totalReactionRate = totalSynthesisRate + conversionRate + fragmentationRate

    // Sampling reaction interval and reaction selection random value
    val timeInterval = ExponentialDistribution(1 / totalReactionRate).sample()
    this.elapsedTime += timeInterval
    val reactionSelector = UniformRealDistribution(0.0, 1.0).sample() * totalReactionRate

    // Reaction definition
    if (reactionSelector < synthesisSup35) {
        // Synthesis of Sup35
        this.sup35MoleculeCount += 1
    } else if (reactionSelector < totalSynthesisRate) {
        // Synthesis of Hsp104
        this.hsp104MoleculeCount += 1
    } else if (reactionSelector < (totalSynthesisRate + fragmentationRate)) {
        // Aggregate fragmentation with a separate function
        val cutSite = (1..fragmentSites).random()
        val finalAggregates = mutableListOf<Int>()

        var totalSites = 0
        var cutMade = false
        for (prionAggregate in this.aggregateSizes) {
            if (totalSites + prionAggregate - 1 < cutSite || cutMade) {
                finalAggregates.add(prionAggregate)
                totalSites += prionAggregate - 1
            } else {
                val leftSize = (cutSite - totalSites)
                val rightSize = prionAggregate - leftSize
                for (aggPiece in listOf<Int>(leftSize, rightSize)) {
                    if (aggPiece < minAggSize) {
                        this.sup35MoleculeCount += aggPiece
                    } else {
                        finalAggregates.add(aggPiece)
                    }
                }
                cutMade = true
            }
        }
        this.aggregateSizes = finalAggregates
        if (this.aggregateSizes.size == 0) this.hasPSI = false
    } else {
        // Conversion of one monomer onto a random aggregate
        val extendingAggregate = (0 until this.aggregateSizes.size).random()
        this.aggregateSizes[extendingAggregate] += 1
        this.sup35MoleculeCount -= 1
    }

}

fun YeastCell.runSimulation(simTime: Double,
                            conversionBeta: Double,
                            fragmentationGamma: Double): MutableList<YeastCell> {

    val offspring = mutableListOf<YeastCell>()
    while (this.elapsedTime < simTime) {
        // println(this.elapsedTime)
        if (this.elapsedTime >= this.nextDivisionTime) {
            // println("Division")
            offspring.add(this.divide())
        } else if (this.hasPSI) {
            // println("Reaction")
            this.nextReaction(conversionBeta = conversionBeta, fragmentationGamma = fragmentationGamma)
        } else {
            // println("Protein levels update")
            this.updateProteinLevel(this.nextDivisionTime)
        }
    }
    return offspring
}
