import java.io.File

// Yeast prion propagation with size threshold and Gillespie's algorithm
// From Derdowski et al., 2012 Science

fun <T, U> cartesianProduct(c1: Collection<T>, c2: Collection<U>): List<Pair<T, U>> {
    return c1.flatMap { lhsElem -> c2.map { rhsElem -> lhsElem to rhsElem } }
}

fun Double.interpolate(endRange: Double, outN: Int): MutableList<Double> {
    val stepSize = (endRange - this) / outN
    var seqElem = this
    val outSeq = mutableListOf<Double>()
    while (seqElem <= endRange) {
        outSeq.add(seqElem)
        seqElem += stepSize
    }
    return outSeq
}

fun main() {

    val betaList = 0.0000065.interpolate(endRange = 0.0003, outN = 16)
    val gammaList = 0.006.interpolate(endRange = 0.013, outN = 25)
    val paramPairs = cartesianProduct(betaList, gammaList)

    val outputFile = File("cultureStats.txt")
    if (outputFile.exists()) { outputFile.delete() }

    paramPairs.parmap { thisPair ->
        val cellCulture = mutableListOf<YeastCell>(YeastCell())
        var cellIndex = 0
        while (cellIndex < cellCulture.size) {
            println("Processing cell $cellIndex ...")
            val newCells = cellCulture[cellIndex].runSimulation(simTime = 1000.0,
                conversionBeta = thisPair.first, fragmentationGamma = thisPair.second)
            cellCulture += newCells
            cellIndex += 1
        }

        val solubleSup35Fraction = cellCulture.map {
            it.sup35MoleculeCount.toDouble() / (it.sup35MoleculeCount + it.aggregateSizes.sum())
        }.average()
        val psiFraction = cellCulture.map { it.hasPSI }.count { true } / cellCulture.size

        outputFile.appendText("${thisPair.first}\t${thisPair.second}\t$psiFraction\t$solubleSup35Fraction\n")

    }
}
