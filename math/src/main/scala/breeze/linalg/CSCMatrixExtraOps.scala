package breeze.linalg

import breeze.linalg.operators.{BinaryOp, OpMulScalar}

object CSCMatrixExtraOps {
  /*
    How to extract data from a CSCMatrix:

    m.colPtrs.sliding(2,1).toList.zipWithIndex.flatMap { case (x, c) =>
      for(ri <- x.head until x.last) yield (c, m.rowIndices(ri), m.data(ri))
    }

   */

  abstract class CSCMatrixCanMulM_M[@specialized (Int, Float, Long, Double) A]
    extends BinaryOp[CSCMatrix[A], CSCMatrix[A], OpMulScalar, CSCMatrix[A]] {

    protected def times(a: A, b: A): A

    protected def zeros  (rows: Int, cols: Int): CSCMatrix[A]
    protected def builder(rows: Int, cols: Int, sz: Int): CSCMatrix.Builder[A]

    final def apply(a: CSCMatrix[A], b: CSCMatrix[A]): CSCMatrix[A] = {
      val rows  = a.rows
      val cols  = a.cols
      require(rows == b.rows, "Matrices must have same number of rows!")
      require(cols == b.cols, "Matrices must have same number of cols!")

      if (cols == 0) return zeros(rows, cols)

      val res     = builder(rows, cols, math.min(a.activeSize, b.activeSize))
      var ci      = 0
      var acpStop = a.colPtrs(0)
      var bcpStop = b.colPtrs(0)
      while (ci < cols) {
        val ci1 = ci + 1
        var acp = acpStop
        var bcp = bcpStop
        acpStop = a.colPtrs(ci1)
        bcpStop = b.colPtrs(ci1)
        while (acp < acpStop && bcp < bcpStop) {
          val ari = a.rowIndices(acp)
          val bri = b.rowIndices(bcp)
          if (ari == bri) {
            val v = times(a.data(acp), b.data(bcp))
            res.add(ari, ci, v)
            acp += 1
            bcp += 1
          } else if (ari < bri) {
            acp += 1
          } else /* ari > bri */ {
            bcp += 1
          }
        }
        ci = ci1
      }

      res.result()
    }
  }

  implicit object CSCMatrixCanMulM_M_Int extends CSCMatrixCanMulM_M[Int] {
    protected def times(a: Int, b: Int) = a * b
    protected def zeros(rows: Int, cols: Int) = CSCMatrix.zeros(rows, cols)
    protected def builder(rows: Int, cols: Int, sz: Int) = new CSCMatrix.Builder(rows, cols, sz)
  }

  implicit object CSCMatrixCanMulM_M_Double extends CSCMatrixCanMulM_M[Double] {
    protected def times(a: Double, b: Double) = a * b
    protected def zeros(rows: Int, cols: Int) = CSCMatrix.zeros(rows, cols)
    protected def builder(rows: Int, cols: Int, sz: Int) = new CSCMatrix.Builder(rows, cols, sz)
  }
}