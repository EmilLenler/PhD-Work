namespace ClassLibrary2;
using static  System.Math;
public class Trap
{
    public double RF_Frequency, kappa, z0, r0;

    //Generator
    public Trap(double RF_Frequency, double kappa, double z0, double r0)
    {
        this.RF_Frequency = RF_Frequency;
        this.kappa = kappa;
        this.r0 = r0;
        this.z0 = z0;
    }
    public void Display()
    {
        Console.WriteLine("This is a trap object");
    }
    public double[] Secular_frequenices(double a, double q)
    {
        double omega_z = RF_Frequency*Sqrt(-a*0.5);
        double omega_r = 0.5*RF_Frequency*Sqrt(0.5*q*q+a);
        double[] freq_array = new double[2];
        freq_array[0] = omega_z;
        freq_array[1] = omega_r;
        return freq_array;
    }
}
