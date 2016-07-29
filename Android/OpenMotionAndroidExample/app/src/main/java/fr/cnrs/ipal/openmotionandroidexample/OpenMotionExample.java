package fr.cnrs.ipal.openmotionandroidexample;


import android.os.Handler;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.widget.TextView;
import java.text.DecimalFormat;
import com.androidplot.xy.SimpleXYSeries;
import com.androidplot.xy.*;
import java.util.Arrays;
import android.graphics.Color;


public class OpenMotionExample extends AppCompatActivity {


    //constants
    private static final int NB_VALUES = 50;

    // Plots
    private XYPlot plot;


    // Arrays of Euler angles values
    private Number[] pitch_values;
    private Number[] roll_values;
    private Number[] yaw_values;

    // SimpleXY Series
    private SimpleXYSeries pitch_series;
    private SimpleXYSeries roll_series;
    private SimpleXYSeries yaw_series;

    // Line Format
    private LineAndPointFormatter format_pitch;
    private LineAndPointFormatter format_roll;
    private LineAndPointFormatter format_yaw;

    // Sensor listener
    private IMUListener mIMUListener;

    // runnable for sensor listener
    private Handler mHandler;
    private final Runnable runnable = new Runnable() {
        public void run() {
            update();
        }
    };


    @Override
    protected void onCreate(Bundle savedInstanceState) {

        super.onCreate(savedInstanceState);

        setContentView(R.layout.activity_open_motion_example);

        // initialization of the Handler
        mHandler = new Handler();

        // initialization of the IMU listener
        mIMUListener = new IMUListener(this.getApplicationContext());

        // initialize our XYPlot reference:
        plot = (XYPlot) findViewById(R.id.plot);

        // create a couple arrays of y-values to plot:
        pitch_values = new Number[]{0.0,0.0,0.0,0.0};
        roll_values = new Number[]{0.0,0.0,0.0,0.0};
        yaw_values = new Number[]{0.0,0.0,0.0,0.0};

        // turn the above arrays into XYSeries':
        pitch_series = new SimpleXYSeries(Arrays.asList(pitch_values),
                SimpleXYSeries.ArrayFormat.Y_VALS_ONLY, "Pitch");

        roll_series = new SimpleXYSeries(Arrays.asList(roll_values),
                SimpleXYSeries.ArrayFormat.Y_VALS_ONLY, "Roll");

        yaw_series = new SimpleXYSeries(Arrays.asList(yaw_values),
                SimpleXYSeries.ArrayFormat.Y_VALS_ONLY, "Yaw");

        // create formatters to use for drawing a series using LineAndPointRenderer
        // and configure them from xml:
        format_pitch = new LineAndPointFormatter();
        format_pitch.setPointLabelFormatter(new PointLabelFormatter(Color.TRANSPARENT));
        format_pitch.configure(getApplicationContext(),R.xml.line_point_formatter_with_labels);
        format_pitch.setInterpolationParams(
                new CatmullRomInterpolator.Params(10, CatmullRomInterpolator.Type.Centripetal));

        format_roll = new LineAndPointFormatter();
        format_roll.setPointLabelFormatter(new PointLabelFormatter(Color.TRANSPARENT));
        format_roll.configure(getApplicationContext(),R.xml.line_point_formatter_with_labels_2);
        format_roll.setInterpolationParams(
                new CatmullRomInterpolator.Params(10, CatmullRomInterpolator.Type.Centripetal));

        format_yaw = new LineAndPointFormatter();
        format_yaw.setPointLabelFormatter(new PointLabelFormatter(Color.TRANSPARENT));
        format_yaw.configure(getApplicationContext(),R.xml.line_point_formatter_with_label_3);
        format_yaw.setInterpolationParams(
                new CatmullRomInterpolator.Params(10, CatmullRomInterpolator.Type.Centripetal));


        // add a new series' to the xyplot:
        plot.addSeries(pitch_series, format_pitch);
        plot.addSeries(roll_series, format_roll);
        plot.addSeries(yaw_series, format_yaw);
        plot.setRangeBoundaries(-180.0,180.0,BoundaryMode.FIXED);


        // reduce the number of range labels
        plot.setTicksPerRangeLabel(3);


        // rotate domain labels 45 degrees to make them more compact horizontally:
        plot.getGraphWidget().setDomainLabelOrientation(-45);


    }


    ///////////////////////////
    // Update
    ///////////////////////////



    private void update(){

        // get IMU output
        mIMUListener.onResume();

        // update label text
        updateText();

        // get Euler angle from OpenMotion
        float[] euler_angle = mIMUListener.getEulerAngle();

        // add last values to the series
        pitch_series.addLast(null,new Double((double)euler_angle[0]));
        roll_series.addLast(null,new Double((double)euler_angle[1]));
        yaw_series.addLast(null,new Double((double)euler_angle[2]));

        // remove eldest values
        if(pitch_series.size() >= NB_VALUES){
            pitch_series.removeFirst();
            roll_series.removeFirst();
            yaw_series.removeFirst();
        }

        // redraw the plot
        plot.redraw();


        mHandler.postDelayed(runnable,100);

    }

    @Override
    public void onResume() {
        super.onResume();
        update();
    }


    @Override
    public void onPause() {
        super.onPause();
        mIMUListener.onPause();
    }

    @Override
    public void onStop() {
        super.onStop();


    }

    @Override
    public void onDestroy() {
        super.onDestroy();
        mIMUListener.onDestroy();
    }


    ///////////////////////////
    //  IMU DISPLAYER
    ///////////////////////////

    public void updateText(){

        DecimalFormat df = new DecimalFormat("0.000##");

        float[] gyro = mIMUListener.getGyroData();
        float[] acc = mIMUListener.getAccData();
        float[] mag = mIMUListener.getMagData();
        float[] quat_openmotion = mIMUListener.getOpenMotionQuat();


        String s_gyro = "Gyroscope output:     [ "+df.format(gyro[0])+" , "+df.format(gyro[1])+" , "+df.format(gyro[2])+" ] ";
        String s_acc = "Accelerometer data: [ "+df.format(acc[0])+" , "+df.format(acc[1])+" , "+df.format(acc[2])+" ] ";
        String s_mag = "Magnetometer data: [ "+df.format(mag[0])+" , "+df.format(mag[1])+" , "+df.format(mag[2])+" ] ";


        String s_quat_om = "Quaternion OpenMotion:  [ "+df.format(quat_openmotion[3])+" , "+df.format(quat_openmotion[0])+" , "+df.format(quat_openmotion[1])+" , "+df.format(quat_openmotion[2])+" ] ";

        TextView t_gyro = (TextView)findViewById(R.id.label_gyro);
        TextView t_acc = (TextView)findViewById(R.id.label_acc);
        TextView t_mag = (TextView)findViewById(R.id.label_mag);
        TextView t_quat_om = (TextView)findViewById(R.id.label_quat_openmotion);

        t_gyro.setText(s_gyro);
        t_acc.setText(s_acc);
        t_mag.setText(s_mag);

        t_quat_om.setText(s_quat_om);

    }



}
