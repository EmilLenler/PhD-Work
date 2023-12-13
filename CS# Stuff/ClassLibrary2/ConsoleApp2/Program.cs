// See https://aka.ms/new-console-template for more information
using ClassLibrary2;
using static System.Math;
public class Runner{
    public static void Main(){
        Trap trap1 = new Trap(5.2*1e6,0.248,3*1e-3,3.5*1e-3);
        trap1.Display();
        double[] sec_freqs = trap1.Secular_frequenices(-0.01,0.2);
        Console.WriteLine((sec_freqs[0]*1e-3/(2*PI)).ToString());
        Console.WriteLine((sec_freqs[1]*1e-3/(2*PI)).ToString());
    }
}