using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Matrix
{
  class MyMatrix
  {

    public static double[][] MatrixCreate(int rows, int cols)
    {
      // Создаем матрицу, полностью инициализированную
      // значениями 0.0. Проверка входных параметров опущена.
      double[][] result = new double[rows][];
      for (int i = 0; i < rows; ++i)
        result[i] = new double[cols]; // автоинициализация в 0.0
      return result;
    }

    public static string MatrixAsString(double[][] matrix)
    {
      string s = "";
      for (int i = 0; i < matrix.Length; ++i)
      {
        for (int j = 0; j < matrix[i].Length; ++j)
          s += matrix[i][j].ToString("F3").PadLeft(8) + " ";
        s += Environment.NewLine;
      }
      return s;
    }

    public static double[][] MatrixProduct(double[][] matrixA, double[][] matrixB)
    {
      // Проверка ошибок, вычисление aRows, aCols, bCols
      int aRows = matrixA.Length; int aCols = matrixA[0].Length;
      int bRows = matrixB.Length; int bCols = matrixB[0].Length;
      if (aCols != bRows)
        throw new Exception("Non-conformable matrices in MatrixProduct");

      double[][] result = MatrixCreate(aRows, bCols);
      Parallel.For(0, aRows, i =>
      {
        for (int j = 0; j < bCols; ++j)
          for (int k = 0; k < aCols; ++k)
            result[i][j] += matrixA[i][k] * matrixB[k][j];
      }
      );
      return result;
    }

    public static double[][] MatrixRandom(int rows, int cols, double minVal, double maxVal, int seed)
    {
      // Возвращаем матрицу со значениями
      // в диапазоне от minVal до maxVal
      Random ran = new Random(seed);
      double[][] result = MatrixCreate(rows, cols);
      for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
          result[i][j] = (maxVal - minVal) * ran.NextDouble() + minVal;
      return result;
    }

    public static double[][] MatrixIdentity(int n)
    {
      double[][] result = MatrixCreate(n, n);
      for (int i = 0; i < n; ++i)
        result[i][i] = 1.0;
      return result;
    }

    public static bool MatrixAreEqual(double[][] matrixA, double[][] matrixB, double epsilon)
    {
      // True, если все значения в A равны
      // соответствующим значениям в B
      int aRows = matrixA.Length;
      int bCols = matrixB[0].Length;
      for (int i = 0; i < aRows; ++i) // каждая строка A и B
        for (int j = 0; j < bCols; ++j) // каждый столбец A и B
          if (Math.Abs(matrixA[i][j] - matrixB[i][j]) > epsilon)
            return false;
      return true;
    }

    public static double[][] MatrixDuplicate(double[][] matrix)
    {
      // Предполагается, что матрица не нулевая
      double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
      for (int i = 0; i < matrix.Length; ++i) // Копирование значений
        for (int j = 0; j < matrix[i].Length; ++j)
          result[i][j] = matrix[i][j];
      return result;
    }

    private static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
    {
      // Разложение LUP Дулитла. Предполагается,
      // что матрица квадратная.
      int n = matrix.Length; // для удобства
      double[][] result = MatrixDuplicate(matrix);
      perm = new int[n];
      for (int i = 0; i < n; ++i) { perm[i] = i; }
      toggle = 1;
      for (int j = 0; j < n - 1; ++j) // каждый столбец
      {
        double colMax = Math.Abs(result[j][j]); // Наибольшее значение в столбце j
        int pRow = j;
        for (int i = j + 1; i < n; ++i)
        {
          if (result[i][j] > colMax)
          {
            colMax = result[i][j];
            pRow = i;
          }
        }
        if (pRow != j) // перестановка строк
        {
          double[] rowPtr = result[pRow];
          result[pRow] = result[j];
          result[j] = rowPtr;
          int tmp = perm[pRow]; // Меняем информацию о перестановке
          perm[pRow] = perm[j];
          perm[j] = tmp;
          toggle = -toggle; // переключатель перестановки строк
        }
        if (Math.Abs(result[j][j]) < 1.0E-20)
          return null;
        for (int i = j + 1; i < n; ++i)
        {
          result[i][j] /= result[j][j];
          for (int k = j + 1; k < n; ++k)
            result[i][k] -= result[i][j] * result[j][k];
        }
      } // основной цикл по столбцу j
      return result;
    }

    private static double[] HelperSolve(double[][] luMatrix, double[] b)
    {
      // Решаем luMatrix * x = b
      int n = luMatrix.Length;
      double[] x = new double[n];
      b.CopyTo(x, 0);
      for (int i = 1; i < n; ++i)
      {
        double sum = x[i];
        for (int j = 0; j < i; ++j)
          sum -= luMatrix[i][j] * x[j];
        x[i] = sum;
      }
      x[n - 1] /= luMatrix[n - 1][n - 1];
      for (int i = n - 2; i >= 0; --i)
      {
        double sum = x[i];
        for (int j = i + 1; j < n; ++j)
          sum -= luMatrix[i][j] * x[j];
        x[i] = sum / luMatrix[i][i];
      }
      return x;
    }

    public static double[][] MatrixInverse(double[][] matrix)
    {
      int n = matrix.Length;
      double[][] result = MatrixDuplicate(matrix);
      int[] perm;
      int toggle;
      double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
      if (lum == null)
        throw new Exception("Unable to compute inverse");
      double[] b = new double[n];
      for (int i = 0; i < n; ++i)
      {
        for (int j = 0; j < n; ++j)
        {
          if (i == perm[j])
            b[j] = 1.0;
          else
            b[j] = 0.0;
        }
        double[] x = HelperSolve(lum, b);
        for (int j = 0; j < n; ++j)
          result[j][i] = x[j];
      }
      return result;
    }

    public static double MatrixDeterminant(double[][] matrix)
    {
      int[] perm;
      int toggle;
      double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
      if (lum == null)
        throw new Exception("Unable to compute MatrixDeterminant");
      double result = toggle;
      for (int i = 0; i < lum.Length; ++i)
        result *= lum[i][i];
      return result;
    }

  }
}
