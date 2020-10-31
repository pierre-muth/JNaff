package jnaff;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class JNaff {

	static int NAFFPoints;
	static double NAFFdt;
	static double[] NAFFData;

	private static int simpleFFT(double[] magnitude2, double[] data, int points) {
		int sizeLimit = points+2;
		int FFTFreqs, i;

		if (sizeLimit < (points+2) ) 
			return 0;

		FastFourierTransformer fastFourierTransformer = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] complex = fastFourierTransformer.transform(data, TransformType.FORWARD);

		FFTFreqs = points/2+1;

		for (i=0; i<FFTFreqs; i++) {
			if (i!=0 && !(i==(FFTFreqs-1) && points%2==0)) {
				magnitude2[i] = Math.pow(complex[i].getReal()*2, 2) + Math.pow(complex[i].getImaginary()*2, 2);
			} else {
				magnitude2[i] = Math.pow(complex[i].getReal(), 2) + Math.pow(complex[i].getImaginary(), 2);
			}

		}
		return FFTFreqs;
	}

	private static double naffFunc(double omega) {
		int i;
		double sum1 = 0, sum2 = 0, cosine, sine;

		for (i=0; i<NAFFPoints; i++) {
			cosine = Math.cos(omega*i*NAFFdt);
			sine   = Math.sin(omega*i*NAFFdt);
			sum1  += cosine*NAFFData[i];
			sum2  += sine*NAFFData[i];
		}

		return Math.pow(sum1, 2) + Math.pow(sum2, 2);
	}

	private static double arithmeticAverage(double[] y, int n) {
		int i;
		double sum = 0;

		for (i=0; i<n; i++) 
			sum += y[i];

		return(sum/n);
	}

	private static int oneDParabolicOptimization(
			int index,
			double[] yReturn, 
			double[] xGuess, 
			double dx,
			double xLower, 
			double xUpper, 
			//			double (*func)(double x, long *invalid),
			double freqCycleLimit, 
			double dxLimit,
			double tolerance, 
			boolean maximize) {

		double maxFactor, fBest, xBest;
		long cycle;
		double x0, x1=0, x2, f0, f1=0, f2;
		double tmp;

		maxFactor = maximize?-1:1;

		x0 = xGuess[index];
		f0 = maxFactor*naffFunc(x0);
		xBest = x0;
		fBest = f0;

		yReturn[index] = maxFactor*f0;

		/* find direction in which function decreases */
		for (cycle=0; cycle<2*freqCycleLimit; cycle++) {
			x1 = x0 + dx;

			if (x1==x0)
				break;
			if (x1>xUpper || x1<xLower)
				return -2;

			f1 = maxFactor*naffFunc(x1);

			if (f1<fBest) {
				fBest = f1;
				xBest = x1;
			}

			if (f1<f0) {
				break;
			}

			dx = dx*(cycle%2==0 ? -1 : -0.5);
		}

		if (x1==x0)
			return 1;
		if (cycle==2*freqCycleLimit) {
			if (Math.abs(dx)<dxLimit) 
				return 1;
			return -3;
		}

		/* take steps until passing through minimum */
		while (true) {
			x2 = x1 + dx;
			if (x2>xUpper || x2<xLower)
				return -4;
			f2 = maxFactor*naffFunc(x2);
			if (f2<fBest) {
				fBest = f2;
				xBest = x2;
			}
			if (f2>f1) 
				break;
			if (x1==x2)
				break;
			x0 = x1; f0 = f1;
			x1 = x2; f1 = f2;
		}

		/* arrange in increasing order */
		if (x0>x2) {
			tmp = x0; x0 = x2; x2 = tmp;
			tmp = f0; f0 = f2; f2 = tmp;
		}

		/* now f0 > f1 and f2 > f1 */
		for (cycle=0; cycle<freqCycleLimit; cycle++) {
			double numer, denom, x3, f3, scale;
			int failed;

			if (x2==x0 || (x2-x0)<dxLimit || ( ( (f2<f0) ? f0 : f2) - f1) < tolerance)
				break;

			/* try parabolic interpolation */
			numer = Math.pow(x1-x0, 2) * (f1-f2) - Math.pow(x1-x2, 2) * (f1-f0);
			denom =    (x1-x0)*(f1-f2) -    (x1-x2)*(f1-f0);
			x3 = x1 - numer/denom/2.0;
			failed = 1;
			scale = x2-x0;

			if (!Double.isInfinite(x3) && 
					x0<x3 && 
					x3<x2 && 
					Math.abs(x3-x0) > 1e-6*scale && 
					Math.abs(x3-x1) > 1e-6*scale &&
					Math.abs(x3-x2) > 1e-6*scale) {

				/* evaluate at parabolic interpolation point */
				failed = 0;
				f3 = maxFactor*naffFunc(x3);

				if (f3 < fBest) {
					fBest = f3;
					xBest = x3;
				}

				if (f3<f1) {
					/* replace point 1 */
					f1 = f3;
					x1 = x3;
				} else if (f2>f0 && f3<f2) {
					/* replace point 2 */
					f2 = f3;
					x2 = x3;
					if (x2<x1) {
						tmp = x1; x1 = x2; x2 = tmp;
						tmp = f1; f1 = f2; f2 = tmp;
					}
				} else if (f2<f0 && f3<f0) {
					/* replace point 0 */
					f0 = f3;
					x0 = x3;
					if (x0>x1) {
						tmp = x0; x0 = x1; x1 = tmp;
						tmp = f0; f0 = f1; f1 = tmp;
					}
				} else failed = 1;

			}

			if (failed == 1) {
				long right, other;
				for (other=0; other<2; other++) {
					/* try dividing one of the intervals */
					if (Math.abs(x0-x1) < Math.abs(x1-x2)) {
						if (other==0) {
							x3 = (x1+x2)/2;
							right = 1;
						} else {
							x3 = (x0+x1)/2;
							right = 0;
						}
					} else {
						if (other==0) {
							x3 = (x0+x1)/2;
							right = 0;
						} else {
							x3 = (x1+x2)/2;
							right = 1;
						}
					}
					f3 = maxFactor*naffFunc(x3);
					if (f3<fBest) {
						fBest = f3;
						xBest = x3;
					}
					if (f3<f1) {
						f1 = f3;
						x1 = x3;
						break;
					}
					if (right!=0 && f3<f2) {
						/* replace point 2 */
						f2 = f3;
						x2 = x3;
						if (x2<x1) {
							tmp = x1; x1 = x2; x2 = tmp;
							tmp = f1; f1 = f2; f2 = tmp;
						}
						break;
					} else if (right==0 && f3<f0) {
						/* replace point 0 */
						f0 = f3;
						x0 = x3;
						if (x0>x1) {
							tmp = x0; x0 = x1; x1 = tmp;
							tmp = f0; f0 = f1; f1 = tmp;
						}
						break;
					}
				}
			}
		}

		yReturn[index] = maxFactor*fBest;
		xGuess[index] = xBest;
		return 1;

	}

	static long calculatePhaseAndAmplitudeFromFreq (
			double[] hanning, 
			int points, 
			double NAFFdt, 
			double frequency, 
			double t0,
			int index,
			double[] phase, 
			double[] amplitude, 
			double[] significance, 
			double[] cosine, 
			double[] sine) {

		int i;
		double sum_ee1, sum_ee2, sum_ef1, sum_ef2;
		double sum1=0, sum2=0, freq0;

		freq0 = frequency;

		sum_ee1 = sum_ee2 = sum_ef1 = sum_ef2 = 0;

		for (i=0; i<points; i++) {
			cosine[i] =  Math.cos(freq0*i*NAFFdt);
			sine[i]   =  Math.sin(freq0*i*NAFFdt);
			/* this gives normalization of the overlap sums */
			sum_ee1 += Math.pow(cosine[i], 2)*hanning[i];
			sum_ee2 += Math.pow(sine[i], 2)*hanning[i];
			/* these are the overlap sums */
			sum_ef1 += cosine[i]*NAFFData[i];
			sum_ef2 += sine[i]*NAFFData[i];
		}

		for (i=0; i<points; i++) 
			sum1 += Math.pow(NAFFData[i], 2);

		for (i=0; i<points; i++)
			NAFFData[i] -= (sum_ef1/sum_ee1*cosine[i] + sum_ef2/sum_ee2*sine[i])*hanning[i];

		for (i=0; i<points; i++) 
			sum2 += Math.pow(NAFFData[i], 2);

		if (sum1>0)
			significance[index] = sum2/sum1;
		else
			significance[index] = -1;

		freq0 = frequency / Math.PI*2;

		amplitude[index] = Math.sqrt( Math.pow(sum_ef1/sum_ee1, 2) + Math.pow(sum_ef2/sum_ee2, 2) );

		/* compute the phase and ensure that it is in the range [-PI, PI] */
		phase[index] = ( Math.atan2(-sum_ef2/sum_ee2, sum_ef1/sum_ee1) + freq0*t0*Math.PI*2 % Math.PI*2);

		if (phase[index] < -Math.PI)
			phase[index] += Math.PI*2;

		if (phase[index] > Math.PI)
			phase[index] -= Math.PI*2;

		return 0;
	}


	public static int performNAFF(
			double[] frequency,   /* return or input frequencies */
			double[] amplitude,   /* return amplitudes */
			double[] phase,       /* return phases */
			double[] significance,/* return "significance" (0=very, large=not) */
			double t0,            /* initial "time" or independent value */
			double dt,            /* spacing of time points */
			double[] data,        /* data to be analyzed */
			int points,           /* number of data points */
			/* these control termination of the iteration for frequencies: */
			/* min acceptable contribution of frequency */
			double fracRMSChangeLimit, 
			/* maximum number of frequencies */
			long maxFrequencies,
			/* these control the accuracy of each frequency: */
			/* maximum iteractions of parabolic optimizer */
			double freqCycleLimit, 
			/* acceptable fractional accuracy of frequency */
			double fracFreqAccuracyLimit,
			/* search only for frequencies between these limits */
			double lowerFreqLimit, 
			double upperFreqLimit
			){

		double rmsOrig, rmsLast, rmsNow, mean;
		int i, freqsFound = 0, FFTFreqs;
		double wStart, freqSpacing;
		int iBest, code, trys;
		double scale, maxMag2;
		double[] sine = new double[points];
		double[] cosine = new double[points];
		double[] magnitude2 = new double[points];
		double[] hanning = new double[points];
		NAFFData = new double[points];

		if (points<2) {
			return -1;
		}

		freqSpacing = 1./(points*dt);
		NAFFdt = dt;

		/* subtract off mean and apply the Hanning window */
		mean = arithmeticAverage(data, points);
		for (i=0; i<points; i++) {
			hanning[i]  = (1 - Math.cos(Math.PI*2*i/(points-1.0)))/2;
			NAFFData[i] = (data[i]-mean)*hanning[i];
		}

		rmsOrig = 0;
		if (fracRMSChangeLimit != 0) {
			for (i=0; i<points; i++)
				rmsOrig += Math.pow(NAFFData[i], 2);
			rmsOrig = Math.sqrt(rmsOrig/points);
		}
		rmsLast = rmsOrig;

		FFTFreqs = points/2-1;
		NAFFPoints = points;
		for (i=0; i<maxFrequencies; i++) 
			amplitude[i] = phase[i] = significance[i] = -1;

		for (i=0; i<maxFrequencies; i++) 
			frequency[i] = -1;

		while (freqsFound < maxFrequencies) {
			simpleFFT(magnitude2, NAFFData, points);
			maxMag2 = 0;
			iBest = 0;
			for (i=0; i < FFTFreqs; i++) {
				if (magnitude2[i] > maxMag2) {
					if (lowerFreqLimit < upperFreqLimit && 
							(i*freqSpacing<lowerFreqLimit || i*freqSpacing>upperFreqLimit))
						continue;
					iBest = i;
					maxMag2 = magnitude2[i];
				}
			}
			
			if (iBest==0)
				break;
			
			wStart = frequency[freqsFound] = iBest*freqSpacing*Math.PI*2;
			scale = naffFunc(wStart);
			for (trys = 0; trys < 2; trys++) {
				code = oneDParabolicOptimization(
						freqsFound,
						amplitude, 
						frequency,
						Math.PI*2*freqSpacing, 
						0.0, 
						Math.PI/dt, 
						freqCycleLimit, 
						fracFreqAccuracyLimit*Math.PI/dt, 
						0.0, 
						true);
				if (code<0) {
					/* amplitude[freqsFound] = frequency[freqsFound] = -1; */
					frequency[freqsFound] = wStart;
					amplitude[freqsFound] = maxMag2;
					break;
				}
			}
			
			calculatePhaseAndAmplitudeFromFreq(
					hanning, 
					points, 
					NAFFdt, 
					frequency[freqsFound], 
					t0,
					freqsFound,
					phase, 
					amplitude,
					significance, 
					cosine, 
					sine);
			
			frequency[freqsFound] /= Math.PI*2;
			freqsFound ++;

			rmsNow = 0;
			
			if (fracRMSChangeLimit != 0) {
				/* determine if residual is too small to bother with */
				for (i=0; i<points; i++) 
					rmsNow += Math.pow(NAFFData[i], 2);
				
				rmsNow = Math.sqrt(rmsNow/points);
				
				if ((rmsLast-rmsNow)/rmsOrig < fracRMSChangeLimit)
					break;
			}
			
			rmsLast = rmsNow;
			
		}
		return freqsFound;
		
	}
	
	public static void main(String[] args) {
		int maxFrequencies = 8;
		double[] frequency = new double[maxFrequencies];
		double[] amplitude = new double[maxFrequencies];
		double[] phase = new double[maxFrequencies];
		double[] significance = new double[maxFrequencies];
		double t0 = 0.0;
		double dt = 1.0;
		int points = 256;
		double[] data = new double[points];        /* data to be analyzed */
		/* these control termination of the iteration for frequencies: */
		/* min acceptable contribution of frequency */
		double fracRMSChangeLimit = 0.0; 
		/* maximum number of frequencies */
		/* these control the accuracy of each frequency: */
		/* maximum iteractions of parabolic optimizer */
		double freqCycleLimit = 100; 
		/* acceptable fractional accuracy of frequency */
		double fracFreqAccuracyLimit = 0.00001;
		/* search only for frequencies between these limits */
		double lowerFreqLimit = 0; 
		double upperFreqLimit = points;
		
		for (int i = 0; i < data.length; i++) {
			data[i] = (Math.sin( ((Math.PI*2)/points) *4*i ) /4 )+
					(Math.cos( ((Math.PI*2)/points) *7.77*i )/8 ) +
					((Math.random()-.5) ) ;
			System.out.println(i+", "+data[i]);
		}
		
		performNAFF(
				frequency, 
				amplitude, 
				phase, 
				significance, 
				t0, 
				dt, 
				data, 
				points, 
				fracRMSChangeLimit, 
				maxFrequencies, 
				freqCycleLimit, 
				fracFreqAccuracyLimit, 
				lowerFreqLimit, 
				upperFreqLimit);
		
		for (int i = 0; i < maxFrequencies; i++) {
			System.out.println(i+", "+ frequency[i]*points +", "+ amplitude[i] +", "+ phase[i] +", "+ significance[i] );
			
		}
		
		
		
	}


}
