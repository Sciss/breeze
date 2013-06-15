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

    protected def zeros  (rows: Int, cols: Int         ): CSCMatrix        [A]
    protected def builder(rows: Int, cols: Int, sz: Int): CSCMatrix.Builder[A]

    final def apply(a: CSCMatrix[A], b: CSCMatrix[A]): CSCMatrix[A] = {
      val rows  = a.rows
      val cols  = a.cols
      require(rows == b.rows, "Matrices must have same number of rows!")
      require(cols == b.cols, "Matrices must have same number of cols!")

      if (cols == 0) return zeros(rows, cols)

      val res     = builder(rows, cols, math.min(a.activeSize, b.activeSize))
      var ci      = 0             // column index [0 ... cols)
      var apStop  = a.colPtrs(0)  // pointer into row indices and data
      var bpStop  = b.colPtrs(0)  // pointer into row indices and data
      while (ci < cols) {
        val ci1 = ci + 1
        var ap  = apStop
        var bp  = bpStop
        apStop = a.colPtrs(ci1)
        bpStop = b.colPtrs(ci1)
        while (ap < apStop && bp < bpStop) {
          val ari = a.rowIndices(ap)  // row index [0 ... rows)
          val bri = b.rowIndices(bp)
          if (ari == bri) {           // column and row match, this cell goes into result matrix
            val v = times(a.data(ap), b.data(bp))
            res.add(ari, ci, v)
            ap += 1
            bp += 1
          } else if (ari < bri) {     // next b row starts further down, therefore increase a pointer
            ap += 1
          } else /* ari > bri */ {    // next a row starts further down, therefore increase b pointer
            bp += 1
          }
        }
        ci = ci1
      }

      res.result()
    }
  }

  implicit object CSCMatrixCanMulM_M_Int extends CSCMatrixCanMulM_M[Int] {
    protected def times(a: Int, b: Int) = a * b

    protected def zeros  (rows: Int, cols: Int         ) =     CSCMatrix.zeros  (rows, cols    )
    protected def builder(rows: Int, cols: Int, sz: Int) = new CSCMatrix.Builder(rows, cols, sz)
  }

  implicit object CSCMatrixCanMulM_M_Double extends CSCMatrixCanMulM_M[Double] {
    protected def times(a: Double, b: Double) = a * b

    protected def zeros  (rows: Int, cols: Int         ) =     CSCMatrix.zeros  (rows, cols    )
    protected def builder(rows: Int, cols: Int, sz: Int) = new CSCMatrix.Builder(rows, cols, sz)
  }
}