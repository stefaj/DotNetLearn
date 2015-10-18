using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DotNetLearn.Mathematics
{
    public static class Matrix
    {
        public static double[] Log(double[] mat)
        {
            double[] newMat = new double[mat.Length];
            for (int i = 0; i < mat.Length; i++)
                newMat[i] = System.Math.Log(mat[i]);
            return newMat;
        }
        public static double[,] Log(double[,] mat)
        {
            double[,] newMat = new double[mat.GetLength(0),mat.GetLength(1)];
            for (int i = 0; i < mat.GetLength(0); i++)
                for (int j = 0; j < mat.GetLength(1); j++)
                    newMat[i, j] = System.Math.Log(mat[i, j]);
            return newMat;
        }
        public static double[,,] Log(double[, ,] mat)
        {
            double[, ,] newMat = new double[mat.GetLength(0), mat.GetLength(1), mat.GetLength(2)];
            for (int i = 0; i < mat.GetLength(0); i++)
                for (int j = 0; j < mat.GetLength(1); j++)
                    for (int k = 0; k < mat.GetLength(2); k++)
                        newMat[i, j, k] = System.Math.Log(mat[i, j, k]);

            return newMat;
        }

        public static double[] DeLog(double[] mat)
        {
            double[] newMat = new double[mat.Length];
            for (int i = 0; i < mat.Length; i++)
                newMat[i] = System.Math.Exp(mat[i]);
            return newMat;
        }
        public static double[,] DeLog(double[,] mat)
        {
            double[,] newMat = new double[mat.GetLength(0), mat.GetLength(1)];
            for (int i = 0; i < mat.GetLength(0); i++)
                for (int j = 0; j < mat.GetLength(1); j++)
                    newMat[i, j] = System.Math.Exp(mat[i, j]);
            return newMat;
        }
        public static double[,,] DeLog(double[, ,] mat)
        {
            double[,,] newMat = new double[mat.GetLength(0), mat.GetLength(1), mat.GetLength(2)];
            for (int i = 0; i < mat.GetLength(0); i++)
                for (int j = 0; j < mat.GetLength(1); j++)
                    for (int k = 0; k < mat.GetLength(2); k++)
                        newMat[i, j, k] = System.Math.Exp(mat[i, j, k]);
            return newMat;
        }
    }
}
