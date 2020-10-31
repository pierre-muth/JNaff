package tests;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import cern.jdve.Chart;
import cern.jdve.ChartInteractor;
import cern.jdve.data.DefaultDataSet;
import cern.jdve.data.DefaultDataSource;
import cern.jdve.renderer.PolylineChartRenderer;
import cern.jdve.renderer.ScatterChartRenderer;
import jnaff.JNaff;


public class NaffChart extends JPanel {
	private DefaultDataSet magDataSet = new DefaultDataSet("mag", new double[]{0,1}, new double[]{0,1});
	private DefaultDataSet rawDataSet = new DefaultDataSet("data", new double[]{0,1}, new double[]{0,1});
	private DefaultDataSet naffDataSet = new DefaultDataSet("naff found", new double[]{0,1}, new double[]{0,1});
	
	public NaffChart() {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				buildGUI();
			}
		});
		
		
		int maxFrequencies = 32;
		double[] frequency = new double[maxFrequencies];
		double[] amplitude = new double[maxFrequencies];
		double[] phase = new double[maxFrequencies];
		double[] significance = new double[maxFrequencies];
		double t0 = 0.0;
		double dt = 1.0;
		int points = (int) Math.pow(2, 8);
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
		
		double[] xValues = new double[points];
		
		// create data
		for (int i = 0; i < data.length; i++) {
			data[i] = (Math.sin( ((Math.PI*2)/points) *100.555*i )  /4 )
					+((Math.random()-.5) *1.0) ;
			System.out.println(i+", "+data[i]);
		}
		
		// create x values
		for (int i = 0; i < xValues.length; i++) {
			xValues[i] = i;
		}
		
		// FFT of data
		FastFourierTransformer fastFourierTransformer = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] complex = fastFourierTransformer.transform(data, TransformType.FORWARD);
		double[] magnitude = new double[complex.length];
		for (int i = 0; i < complex.length; i++) {
			magnitude[i] = Math.sqrt( complex[i].getReal()*complex[i].getReal() + complex[i].getImaginary()*complex[i].getImaginary() );
		}
		
		// naff of data
		JNaff.performNAFF(
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
			frequency[i]*=points;
			System.out.println(i+", "+ frequency[i] +", "+ amplitude[i] +", "+ phase[i] +", "+ significance[i] );
			
		}
		
		rawDataSet.set(xValues, data);
		magDataSet.set(xValues, magnitude);
		naffDataSet.set(frequency, amplitude);
		
		
	}
	
	private void buildGUI() {
		Font font = new Font(Font.DIALOG, Font.PLAIN, 10);
		setLayout(new BorderLayout());
		setPreferredSize(new Dimension(800, 600));
		
		Chart frequencyChart = new Chart();
		Chart timeChart = new Chart();
		frequencyChart.setPreferredSize(new Dimension(800, 300));
		timeChart.setPreferredSize(new Dimension(800, 300));
		
		frequencyChart.addInteractor(ChartInteractor.ZOOM);
		frequencyChart.addInteractor(ChartInteractor.DATA_PICKER);
		
		PolylineChartRenderer rawRenderer = new PolylineChartRenderer();
		DefaultDataSource rawDataSource = new DefaultDataSource();
		rawDataSource.addDataSet(rawDataSet);
		rawRenderer.setDataSource(rawDataSource);
		
		ScatterChartRenderer naffRenderer = new ScatterChartRenderer();
		DefaultDataSource naffDataSource = new DefaultDataSource();
		naffDataSource.addDataSet(naffDataSet);
		naffRenderer.setDataSource(naffDataSource);
		
		PolylineChartRenderer magRenderer = new PolylineChartRenderer();
		DefaultDataSource magDataSource = new DefaultDataSource();
		magDataSource.addDataSet(magDataSet);
		magRenderer.setDataSource(magDataSource);
		
		timeChart.addRenderer(rawRenderer);
		
		frequencyChart.addYAxis(true, false);
		frequencyChart.addRenderer(0, magRenderer);
		frequencyChart.addRenderer(1, naffRenderer);
		
		add(frequencyChart, BorderLayout.CENTER);
		add(timeChart, BorderLayout.NORTH);
		
	}

	public static void main(String[] args) {
		JFrame frame = new JFrame("Jnaff Demo");
		frame.setLayout(new BorderLayout());
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		NaffChart main = new NaffChart();
		frame.add(main, BorderLayout.CENTER );
		frame.pack();
		frame.setVisible(true);
	}

}
