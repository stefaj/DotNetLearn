using DotNetLearn.Mathematics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


// TODO
// Fix underflows
// Add additional learning technique
// Support multiple observations training

namespace DotNetLearn.Markov
{
    public class HiddenMarkovModel
    {

        //   double[,] transition;
        //    double[,] emission;
        //   double[] pi;


        double[,] logTransition;
        double[,] logEmission;
        double[] logPi;

        //double[,] alpha;
        //double[,] beta;
        //double[, ,] digamma;
        //double[,] gamma;


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
        public static HiddenMarkovModel FromMatrices(double[,] A, double[,] B, double[] pi, bool isLog = false)
        {
            HiddenMarkovModel model = new HiddenMarkovModel(A.GetLength(0), B.GetLength(1));

            if (!isLog)
            {
                model.logTransition = MathE.elnify(A);
                model.logEmission = MathE.elnify(B);
                model.logPi = MathE.elnify(pi);
            }
            else
            {
                model.logTransition = A;
                model.logEmission = B;
                model.logPi = pi;
            }
            return model;
        }

        public HiddenMarkovModel(int states, int observationSymbols)
        {
            this.logTransition = new double[states, states];
            this.logEmission = new double[states, observationSymbols];
            this.logPi = new double[states];
            this.states = states;
            this.observationSymbols = observationSymbols;
            this.random = new DotNetLearn.Statistics.SimpleRNG();
            int N = states;
            int M = observationSymbols;
        }

        public void UniformRandomPriors()
        {

            logPi = MathE.elnify(NormalMass(logPi.Length));

            for (int i = 0; i < logTransition.GetLength(0); i++)
            {
                double[] trans = MathE.elnify(NormalMass(states));
                for (int j = 0; j < logTransition.GetLength(0); j++)
                {
                    logTransition[i, j] = trans[j];
                }
            }

            for (int i = 0; i < logEmission.GetLength(0); i++)
            {
                double[] ems = MathE.elnify(NormalMass(observationSymbols));
                for (int j = 0; j < logEmission.GetLength(1); j++)
                {
                    logEmission[i, j] = ems[j];
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



        public double[,] CalculateLogGamma(int[] O)
        {
            int T = O.Length;
            Initialize(T);
            var alpha = LogSpaceForward(O);
            var beta = LogSpaceBackward(O);
            var gamma = LogSpaceGamma(O, alpha, beta);
            return gamma;

        }


        public double[,] CalculateGamma(int[] O)
        {
            return MathE.eexpify(CalculateLogGamma(O));
        }

        private void Initialize(int T)
        {
            int N = states;
            int M = observationSymbols;





        }


        // ln(a_t(i))
        private double[,] LogSpaceForward(int[] O)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            var logAlpha = new double[T, N];

            #region Alpha Pass

            for (int i = 0; i < N; i++)
            {
                //logAlpha[0, i] = MathE.elnproduct(MathE.eln(pi[i]), MathE.eln(emission[i, O[0]]));
                logAlpha[0, i] = MathE.elnproduct(logPi[i], logEmission[i, O[0]]);
            }

            //compute a_t(i)
            for (int t = 1; t < T; t++)
            {
                for (int i = 0; i < N; i++)
                {
                    double logalpha = MathE.LOGZERO;
                    //alpha[t, i] = 0;
                    for (int j = 0; j < N; j++)
                    {
                        //logalpha = MathE.elnsum(logalpha, MathE.elnproduct(logAlpha[t - 1, j], MathE.eln(transition[j, i])));
                        logalpha = MathE.elnsum(logalpha, MathE.elnproduct(logAlpha[t - 1, j], logTransition[j, i]));
                    }

                    // logAlpha[t, i] = MathE.elnproduct(logalpha, MathE.eln(emission[i, O[t]]));
                    logAlpha[t, i] = MathE.elnproduct(logalpha, logEmission[i, O[t]]);
                }


            }

            return logAlpha;

            #endregion

        }

        // ln(b_t(i))
        private double[,] LogSpaceBackward(int[] O)
        {

            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            var logBeta = new double[T, N];

            #region Beta Pass

            for (int i = 0; i < N; i++)
            {
                logBeta[T - 1, i] = 0;
            }

            // Beta pass
            for (int t = T - 2; t >= 0; t--)
            {
                for (int i = 0; i < N; i++)
                {
                    double logbeta = MathE.LOGZERO;

                    for (int j = 0; j < N; j++)
                    {
                        // logbeta = MathE.elnsum(logbeta, MathE.elnproduct(MathE.eln(transition[i, j]),
                        //    MathE.elnproduct(emission[j, O[t + 1]], logBeta[t + 1, j])));
                        // mISTAKE ??? MathE.elnproduct(emission[j, O[t + 1]],...

                        logbeta = MathE.elnsum(logbeta, MathE.elnproduct(logTransition[i, j],
                            MathE.elnproduct(logEmission[j, O[t + 1]], logBeta[t + 1, j])));


                    }
                    logBeta[t, i] = logbeta;
                }
            }
            #endregion

            return logBeta;
        }

        private double[,] LogSpaceGamma(int[] O, double[,] logAlpha, double[,] logBeta)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            var logGamma = new double[T, N];

            #region Calculating digamma and gamma
            for (int t = 0; t < T; t++)
            {
                double normalizer = MathE.LOGZERO;
                for (int i = 0; i < N; i++)
                {
                    logGamma[t, i] = MathE.elnproduct(logAlpha[t, i], logBeta[t, i]);
                    normalizer = MathE.elnsum(normalizer, logGamma[t, i]);
                }

                for (int i = 0; i < N; i++)
                {
                    logGamma[t, i] = MathE.elnproduct(logGamma[t, i], -normalizer);
                }
            }

            #endregion

            return logGamma;
        }

        private double[, ,] LogSpaceEta(int[] O, double[,] logAlpha, double[,] logBeta)
        {
            int N = states;
            int M = observationSymbols;
            int T = O.Length;

            var logEta = new double[T, N, N];

            #region Calculating digamma and gamma
            for (int t = 0; t < T - 1; t++)
            {
                double normalizer = MathE.LOGZERO;
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {

                        //logEta[t, i, j] = MathE.elnproduct(logAlpha[t, i], MathE.elnproduct(MathE.eln(transition[i, j]),
                        //   MathE.elnproduct(MathE.eln(emission[j, O[t + 1]]), logBeta[t + 1, j])));

                        /*logEta[t, i, j] = MathE.elnproduct(logAlpha[t, i], MathE.elnproduct(logTransition[i, j],
                            MathE.elnproduct(logEmission[j, O[t + 1]], logBeta[t + 1, j])));
                         */
                        logEta[t, i, j] = MathE.elnproduct(logAlpha[t, i], logTransition[i, j], logEmission[j, O[t + 1]], logBeta[t + 1, j]);
                        normalizer = MathE.elnsum(normalizer, logEta[t, i, j]);

                    }
                }

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        logEta[t, i, j] = MathE.elnproduct(logEta[t, i, j], -normalizer);
                    }
                }
            }

            #endregion


            return logEta;
        }

        private void UpdateModel(int[][] Obs)
        {
            int N = states;
            int M = observationSymbols;

            int K = Obs.GetLength(0);





            Dictionary<int, double[,]> alphas = new Dictionary<int, double[,]>();
            Dictionary<int, double[,]> betas = new Dictionary<int, double[,]>();
            Dictionary<int, double[,]> gammas = new Dictionary<int, double[,]>();
            Dictionary<int, double[, ,]> etas = new Dictionary<int, double[, ,]>();

            for (int k = 0; k < K; k++)
            {
                int[] O = Obs[k];
                int T = O.Length;

                var logAlpha = LogSpaceForward(O);
                var logBeta = LogSpaceBackward(O);
                var logGamma = LogSpaceGamma(O, logAlpha, logBeta);
                var logEta = LogSpaceEta(O, logAlpha, logBeta);

                alphas[k] = logAlpha;
                betas[k] = logBeta;
                gammas[k] = logGamma;
                etas[k] = logEta;
            }

            #region Re-estimation of parameters




            // re-sestimate pi



            for (int i = 0; i < N; i++)
            {
                double logpi = MathE.LOGZERO;
                for (int k = 0; k < K; k++)
                {

                    var logGamma = gammas[k];
                    int T = Obs[k].Length;

                    logpi = MathE.elnsum(logpi, logGamma[0, i]);
                }
                logPi[i] = MathE.elnproduct(logpi, -K);

            }



            // re-estimate A
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    double numer = MathE.LOGZERO;
                    double denom = MathE.LOGZERO;

                    for (int k = 0; k < K; k++)
                    {
                        var logEta = etas[k];
                        var logGamma = gammas[k];
                        int T = Obs[k].Length;

                        for (int t = 0; t < T - 1; t++)
                        {
                            numer = MathE.elnsum(numer, logEta[t, i, j]);
                            denom = MathE.elnsum(denom, logGamma[t, i]);
                        }
                    }
                    logTransition[i, j] = MathE.elnproduct(numer, -denom);
                }
            }



            // re-estimate B
            for (int j = 0; j < N; j++)
            {
                for (int ki = 0; ki < M; ki++)
                {
                    double numer = MathE.LOGZERO;
                    double denom = MathE.LOGZERO;
                    for (int k = 0; k < K; k++)
                    {
                        var logEta = etas[k];
                        var logGamma = gammas[k];
                        int T = Obs[k].Length;

                        for (int t = 0; t < T; t++)
                        {
                            if (Obs[k][t] == ki)
                            {
                                numer = MathE.elnsum(numer, logGamma[t, j]);
                            }
                            denom = MathE.elnsum(denom, logGamma[t, j]);
                        }
                    }
                    logEmission[j, ki] = MathE.elnproduct(numer, -denom);
                }
            }

            #endregion
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
            var logAlpha = LogSpaceForward(O);
            var logBeta = LogSpaceBackward(O);
            var logGamma = LogSpaceGamma(O, logAlpha, logBeta);
            var logEta = LogSpaceEta(O, logAlpha, logBeta);

            double[] stateSeq = new double[T];
            for (int t = 0; t < T; t++)
            {
                double max = 0;
                double bestState = 0;
                for (int i = 0; i < N; i++)
                {
                    if (logGamma[t, i] > max)
                    {
                        max = logGamma[t, i];
                        bestState = i;
                    }
                }

                stateSeq[t] = bestState;
            }

            return stateSeq;
        }

        public double EvaluateProb(int[] O)
        {
            var alpha = LogSpaceForward(O);
            return EvaluateProb(O.Length, alpha);
        }

        private double EvaluateProb(int T, double[,] logAlpha)
        {
            int N = states;

            double prob = 0;
            for (int i = 0; i < N; i++)
            {
                prob += MathE.eexp(logAlpha[T - 1, i]);
            }
            return prob;
        }

        // Independence assumption
        public double EvaluateProb(int[][] O)
        {
            int K = O.GetLength(0);
            double prob = 1;

            for (int k = 0; k < K; k++)
            {
                prob *= EvaluateProb(O[k]);
            }

            return prob;
        }


        public void Train(int[][] O)
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

                UpdateModel(O);

                //double logProb = EvaluateLogProb(T);

                double logProb = EvaluateProb(O);

                if (oldLogProb > logProb && iter > MaxIterations / 2)
                    break;
                Console.WriteLine(logProb);

                oldLogProb = logProb;


            }



        }


        public int[] PredictObservationSequence(int T)
        {
            var emission = MathE.eexpify(logEmission);
            var transition = MathE.eexpify(logTransition);
            var pi = MathE.eexpify(logPi);

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
            for (int j = 0; j < logEmission.GetLength(1); j++)
            {
                for (int i = 0; i < logEmission.GetLength(0); i++)
                {
                    Console.Write("{0:.00} ", MathE.eexp(logEmission[i, j]));
                }
                Console.WriteLine();
            }
        }

        public void PrintTransitions()
        {
            for (int j = 0; j < logTransition.GetLength(1); j++)
            {
                for (int i = 0; i < logTransition.GetLength(0); i++)
                {
                    Console.Write("{0:.00} ", MathE.eexp(logTransition[i, j]));
                }
                Console.WriteLine();
            }
        }

        public void PrintPriorPi()
        {

            for (int i = 0; i < logPi.GetLength(0); i++)
            {
                Console.Write("{0:.00} ", MathE.eexp(logPi[i]));
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
