using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Matrix
{
  class Program
  {
    static void Main(string[] args)
    {

      double[][] m = MyMatrix.MatrixCreate(3, 3);

      m[0][0] = 1;
      m[0][1] = 2;
      m[0][2] = 3;

      m[1][0] = 0;
      m[1][1] = 77;
      m[1][2] = 3;

      m[2][0] = 11;
      m[2][1] = 0;
      m[2][2] = 14;

      var inv = MyMatrix.MatrixInverse(m);
      Console.WriteLine(MyMatrix.MatrixAsString(m));
      Console.WriteLine("Inverse = " + MyMatrix.MatrixAsString(inv));

    }
  }
}
