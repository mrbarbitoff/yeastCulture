// Yeast prion propagation with size threshold and Gillespie's algorithm
// From Derdowski et al., 2012 Science


fun main(args: Array<String>) {
    println("Hello World!")
    println(args)
    println(totalSynthesisRate)

    var divCount: Int = 1

    var testCdll = yeastCell()
    testCdll.updateProteinLevel(1000.0)
    var daughterCell = testCdll.divide()
    println(testCdll.aggregateSizes)

    val testList = listOf<Double>(0.01, 0.05)
    println(testList.parmap { it * 6 })
}