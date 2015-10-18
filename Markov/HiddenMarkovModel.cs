using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DotNetLearn.Markov
{
    public class HiddenMarkovModel
    {
        double[,] transition;
        double[,] emission;
        double[] pi;

        double[,] alpha;
        double[,] beta;
        double[, ,] digamma;
        double[,] gamma;
        double[] c;

        public int MaxIterations { get; set; }

        int states;
        int observationSymbols;

        DotNetLearn.Statistics.SimpleRNG random;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="A">Transition matrix</param>
        /// <param name="B">Emision matrix</param>
        /// <returns></returns>
        public static HiddenMarkovModel FromMatrices(double[,] A, double[,] B, double[] pi)
        {
            HiddenMarkovModel model = new HiddenMarkovModel(A.GetLength(0), B.GetLength(1));
            model.transition = A;
            model.emission = B;
            model.pi = pi;
            return model;
        }

        public HiddenMarkovModel(int states, int observationSymbols)
        {
            this.transition = new double[states, states];
            this.emission = new double[states, observationSymbols];
            this.pi = new double[states];
            this.states = states;
            this.observationSymbols = observationSymbols;
            this.random = new DotNetLearn.Statistics.SimpleRNG();
            int N = states;
            int M = observationSymbols;
        }

        public void UniformRandomPriors()
        {

            pi = NormalMass(pi.Length);

            for (int i = 0; i < transition.GetLength(0); i++)
            {
                double[] trans = NormalMass(states);
                for (int j = 0; j < transition.GetLength(0); j++)
                {

                    transition[i, j] = trans[j];
                }
            }

            for (int i = 0; i < emission.GetLength(0); i++)
            {
                double[] ems = NormalMass(observationSymbols);
                for (int j = 0; j < emission.GetLength(1); j++)
                {

                    emission[i, j] = ems[j];
                }
            }
        }

        /// <summary>
        /// Returns a array that is almost uniform and sums to 1
        /// </summary>
        /// <param name="n">size of returned array</param>
        /// <param name="std">Standard deviation</param>
        /// <returns></returns>
        private double[] NormalMass(int n, double std = 0.01)
        {
            double[] seq = new double[n];

            double uni = 1 / (double)n;
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                double r = random.GetNormal(uni, std);
                seq[i] = r;
                sum += r;
            }

            // normalize
            for (int i = 0; i < n; i++)
            {
                seq[i] /= sum;
            }

            return seq;
        }

        public double[,] CalculateGamma(int[] O)
        {
            int T = O.Length;
            Initialize(T);
            ForwardPass(O);
            BackwardPass(O);
            GammaPass(O);
            return gamma;
        }

        private void Initialize(int T)
        {
            int N = states;
            int M = observationSymbols;


            alpha = new double[T, N];
            beta = new double[T, N];
            digamma = new double[T, N, N];
            gamma = new double[T, N];

            c = new double[T];
        }

        private void ForwardPass(int[] O)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            #region Alpha Pass
            //compute alpha a_0(i)
            c[0] = double.Epsilon;
            for (int i = 0; i < N; i++)
            {
                alpha[0, i] = pi[i] * emission[i, O[0]];
                c[0] = c[0] + alpha[0, i];

            }


            //scale  a_0(i)
            c[0] = 1 / c[0];
            for (int i = 0; i < N; i++)
            {
                if (alpha[0, i] != 0)
                    alpha[0, i] = c[0] * alpha[0, i];
            }

            //compute a_t(i)
            for (int t = 1; t < T; t++)
            {
                c[t] = 0;
                for (int i = 0; i < N; i++)
                {

                    alpha[t, i] = 0;
                    for (int j = 0; j < N; j++)
                    {
                        alpha[t, i] = alpha[t, i] + alpha[t - 1, j] * transition[j, i];
                    }

                    alpha[t, i] = alpha[t, i] * emission[i, O[t]];
                    c[t] = c[t] + alpha[t, i];

                }
                //Scale a_t(i)
                c[t] = 1 / c[t];
                for (int i = 0; i < N; i++)
                {
                    double prev_alpha = alpha[t, i];
                    if (alpha[t, i] != 0)
                        alpha[t, i] = c[t] * alpha[t, i];

                }


            }

            #endregion

        }

        private void BackwardPass(int[] O)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            #region Beta Pass
            // Scale B_T-1(i)
            for (int i = 0; i < N; i++)
            {
                beta[T - 1, i] = c[T - 1];
            }

            // Beta pass
            for (int t = T - 2; t >= 0; t--)
            {
                for (int i = 0; i < N; i++)
                {
                    beta[t, i] = 0;
                    for (int j = 0; j < N - 1; j++)
                    {
                        beta[t, i] += transition[i, j] * emission[j, O[t + 1]] * beta[t + 1, j];
                        if (transition[i, j] == 0 || emission[j, O[t + 1]] == 0 || beta[t + 1, j] == 0)
                            beta[t, i] = 0;

                    }
                    //scale B
                    if (beta[t, i] != 0)
                        beta[t, i] = c[t] * beta[t, i];

                }
            }
            #endregion
        }

        private void GammaPass(int[] O)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;


            #region Calculating digamma and gamma
            for (int t = 0; t < T - 1; t++)
            {
                double denom = double.Epsilon;
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        denom += alpha[t, i] * transition[i, j] * emission[j, O[t + 1]] * beta[t + 1, j];
                    }
                }

                for (int i = 0; i < N; i++)
                {
                    gamma[t, i] = 0;
                    for (int j = 0; j < N; j++)
                    {
                        digamma[t, i, j] = (alpha[t, i] * transition[i, j] * emission[j, O[t + 1]] * beta[t + 1, j]) / denom;
                        gamma[t, i] += digamma[t, i, j];
                    }

                }
            }

            // special case for gamma_T-1(i)
            double denom2 = 0;
            for (int i = 0; i < N; i++)
            {
                denom2 += alpha[T - 1, i];
            }
            for (int i = 0; i < N; i++)
            {
                gamma[T - 1, i] = alpha[T - 1, i] / denom2;
            }

            #endregion
        }

        private void ReEstimateModel(int[] O)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            #region Re-estimation of parameters
            //re-estimate pi
            for (int i = 0; i < N; i++)
            {
                pi[i] = gamma[0, i];
            }

            // re-estimate A
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    double numer = 0;
                    double denom = double.Epsilon;
                    for (int t = 0; t < T - 1; t++)
                    {
                        numer += digamma[t, i, j];
                        denom += gamma[t, i];

                    }

                    transition[i, j] = numer / denom;
                }
            }

            // re-estimate B
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    double numer = 0;
                    double denom = double.Epsilon;
                    for (int t = 0; t < T; t++)
                    {
                        if (O[t] == j)
                        {
                            numer += gamma[t, i];
                        }
                        denom += gamma[t, i];
                    }

                    emission[i, j] = numer / denom;

                }
            }

            #endregion
        }

        private double EvaluateLogProb(int T)
        {

            double logProb = 0;
            for (int i = 0; i < T; i++)
            {
                logProb += Math.Log(c[i]);
            }
            return -logProb;
        }


        /// <summary>
        /// Finds the most likely state sequence given an observation
        /// </summary>
        /// <param name="O">Observation</param>
        /// <returns></returns>
        public double[] FindLikelyStateSequence(int[] O)
        {
            int T = O.Length;
            int N = states;

            Initialize(T);
            ForwardPass(O);
            BackwardPass(O);
            GammaPass(O);

            double[] stateSeq = new double[T];
            for (int t = 0; t < T; t++)
            {
                double max = 0;
                double bestState = 0;
                for (int i = 0; i < N; i++)
                {
                    if (gamma[t, i] > max)
                    {
                        max = gamma[t, i];
                        bestState = i;
                    }
                }

                stateSeq[t] = bestState;
            }

            return stateSeq;
        }

        public double EvaluteLogProb(int[] O)
        {
            ForwardPass(O);
            return EvaluateLogProb(O.Length);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="O">Observation Sequence</param>
        public void Train(int[] O)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            #region Initialization

            Initialize(T);


            double oldLogProb = float.NegativeInfinity;
            #endregion

            for (int iter = 0; iter < MaxIterations; iter++)
            {
                ForwardPass(O);
                BackwardPass(O);
                GammaPass(O);
                ReEstimateModel(O);

                double logProb = EvaluateLogProb(T);

                if (oldLogProb > logProb && iter > MaxIterations / 10)
                    break;
                //Console.WriteLine(logProb);

                oldLogProb = logProb;


            }

        }


        public int[] PredictObservationSequence(int T)
        {
            List<int> O = new List<int>();
            Random rand = new Random();

            int[] states_ = new int[states];
            for (int j = 0; j < states; j++)
            {
                states_[j] = j;
            }
            int state = RouletteSelection<int>(states_, pi, rand);


            for (int t = 0; t < T; t++)
            {
                // Get next observation state
                double[] probabilities = new double[observationSymbols];
                int[] observations = new int[observationSymbols];
                bool stop = true;
                for (int j = 0; j < observationSymbols; j++)
                {
                    probabilities[j] = emission[state, j];
                    observations[j] = j;
                    if (probabilities[j] > 0)
                        stop = false;
                }
                if (stop)
                    break;
                int winner = RouletteSelection<int>(observations, probabilities, rand);
                O.Add(winner);

                stop = true;
                // Proceed to next hidden state
                probabilities = new double[states];
                int[] nextStates = new int[states];
                for (int j = 0; j < states; j++)
                {
                    probabilities[j] = transition[state, j];
                    nextStates[j] = j;
                    if (probabilities[j] > 0)
                        stop = false;
                }
                if (stop)
                    break;

                winner = RouletteSelection<int>(nextStates, probabilities, rand);
                state = winner;
            }
            return O.ToArray();
        }

        public void PrintEmission()
        {
            for (int j = 0; j < emission.GetLength(1); j++)
            {
                for (int i = 0; i < emission.GetLength(0); i++)
                {
                    Console.Write("{0:.00} ", emission[i, j]);
                }
                Console.WriteLine();
            }
        }

        public void PrintTransitions()
        {
            for (int j = 0; j < transition.GetLength(1); j++)
            {
                for (int i = 0; i < transition.GetLength(0); i++)
                {
                    Console.Write("{0:.00} ", transition[i, j]);
                }
                Console.WriteLine();
            }
        }

        public void PrintPriorPi()
        {

            for (int i = 0; i < pi.GetLength(0); i++)
            {
                Console.Write("{0:.00} ", pi[i]);
            }
            Console.WriteLine();
        }




        /// <summary>
        /// Returns a random item according to the proportions
        /// </summary>
        /// <typeparam name="B"></typeparam>
        /// <param name="items">Array of values</param>
        /// <param name="proportions">Weight of each item</param>
        /// <param name="rand">Uses own random if none is given</param>
        /// <returns></returns>
        public static B RouletteSelection<B>(B[] items, double[] proportions, Random rand = null)
        {

            if (items.Length != proportions.Length)
                throw new Exception("WTP length mismatch");

            double total_freq = 0;
            foreach (var freq in proportions)
                total_freq += freq;

            // Cumulative probabilities
            double[] probabilities_cum = new double[proportions.Length];
            for (int i = 0; i < items.Length; i++)
            {
                double proportion = proportions[i] / total_freq;

                if (i == 0)
                    probabilities_cum[i] = proportion;
                else
                    probabilities_cum[i] = probabilities_cum[i - 1] + proportion;
            }

            double randomDouble = (double)(rand == null ? StaticRandom.NextDouble() : rand.NextDouble());
            double random_prob = (double)(probabilities_cum[probabilities_cum.Length - 1] * randomDouble);
            int index = Array.BinarySearch<double>(probabilities_cum, random_prob);
            if (index < 0)
                index = -(index + 1);
            if (index == items.Length)
                index = items.Length - 1;
            return items[index];
        }

    }
}
