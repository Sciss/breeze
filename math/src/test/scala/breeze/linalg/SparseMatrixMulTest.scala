package breeze.linalg

object SparseMatrixMulTest extends App {
  import CSCMatrixExtraOps._
  val m1  = CSCMatrix((0, 0, 0), (0, 5, 0), (0, 0, 10), (0, 13, 0))
  println(m1.toDenseMatrix)
  println(s"\ncolPtrs:\n${m1.colPtrs.mkString(", ")}\n\nrowIndices:\n${m1.rowIndices.mkString(", ")}\n")
  val res = m1 :* m1
  println(res.toDenseMatrix)
}