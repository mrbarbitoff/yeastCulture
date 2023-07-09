// Yeast prion propagation with size threshold and Gillespie's algorithm
// From Derdowski et al., 2012 Science


fun main(args: Array<String>) {
    println("Hello World!")
    println(args)

    val testList = listOf<Double>(0.01, 0.05)
    println(testList.parmap { it * 6 })
}