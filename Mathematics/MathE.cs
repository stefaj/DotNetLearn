using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DotNetLearn.Mathematics
{
    public static class MathE
    {
        public const double LOGZERO = double.NaN;

        public static double eexp(double x)
        {
            if (double.IsNaN(x))
                return 0;
            return Math.Exp(x);
        }

        public static double eln(double x)
        {
            if (x == 0)
                return LOGZERO;
            return Math.Log(x);
        }

        public static double elnsum(double eln_x, double eln_y)
        {
            if (double.IsNaN(eln_x))
                return eln_y;
            else if (double.IsNaN(eln_y))
                return eln_x;
            else
            {
                if (eln_x > eln_y)
                {
                    return eln_x + eln(1 + Math.Exp(eln_y - eln_x));
                }
                else
                {
                    return eln_y + eln(1 + Math.Exp(eln_x - eln_y));
                }
            }
        }

        public static double elnproduct(double eln_x, double eln_y)
        {
            if (double.IsNaN(eln_x))
                return LOGZERO;
            if (double.IsNaN(eln_y))
                return LOGZERO;

            return eln_x + eln_y;
        }

        public static double elnproduct(params double[] lns)
        {
            for (int i = 0; i < lns.Length; i++)
            {
                if (double.IsNaN(lns[i]))
                    return LOGZERO;
            }
            double sum = 0;
            for (int i = 0; i < lns.Length; i++)
                sum += lns[i];
            return sum;
        }


        public static double[] eexpify(double[] mat)
        {
            var o = new double[mat.Length];

            for (int i = 0; i < mat.Length; i++)
                o[i] = MathE.eexp(mat[i]);

            return o;
        }

        public static double[,] eexpify(double[,] mat)
        {
            var o = new double[mat.GetLength(0), mat.GetLength(1)];

            for (int i = 0; i < mat.GetLength(0); i++)
                for (int j = 0; j < mat.GetLength(1); j++)
                    o[i, j] = MathE.eexp(mat[i, j]);

            return o;
        }

        public static double[] elnify(double[] mat)
        {
            var o = new double[mat.Length];

            for (int i = 0; i < mat.Length; i++)
                o[i] = MathE.eln(mat[i]);

            return o;
        }

        public static double[,] elnify(double[,] mat)
        {
            var o = new double[mat.GetLength(0), mat.GetLength(1)];

            for (int i = 0; i < mat.GetLength(0); i++)
                for (int j = 0; j < mat.GetLength(1); j++)
                    o[i, j] = MathE.eln(mat[i, j]);

            return o;
        }




    }
}
