param(
    [Parameter(Mandatory=$true)][string]$ImageDirectory,
    [Parameter(Mandatory=$true)][string]$OutputDirectory
)
# Read-only pixel analysis of thesis PDF pages 45-47 rendered at 9000 x 12729.
# Source and axis calibration are documented with the resulting digitized CSVs.
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Drawing
Add-Type -ReferencedAssemblies System.Drawing -TypeDefinition @'
using System;
using System.IO;
using System.Drawing;
using System.Collections.Generic;
using System.Globalization;
public class CocoaFigureDigitizer {
    static bool Dark(Color c) { return c.R<180 && Math.Abs(c.R-c.G)<8 && Math.Abs(c.G-c.B)<8; }
    static bool Red(Color c) { return c.R>180 && c.G<90 && c.B<90; }
    static List<double> Centers(List<int> pixels) {
        var result=new List<double>();
        for(int i=0;i<pixels.Count;) {
            int first=pixels[i], last=first; i++;
            while(i<pixels.Count && pixels[i]<=last+1) { last=pixels[i]; i++; }
            result.Add((first+last)/2.0);
        }
        return result;
    }
    // The central prediction is solid; credible limits are dashed. In a
    // 61-column window retain a connected red component spanning >=90% of
    // the window and crossing its center. Ambiguous/occluded points are omitted.
    static double? SolidY(Bitmap b,int center,int top,int bottom) {
        int left=center-30,width=61,height=bottom-top+1;
        bool[] pixels=new bool[width*height];
        for(int y=0;y<height;y++) for(int x=0;x<width;x++)
            pixels[y*width+x]=Red(b.GetPixel(left+x,top+y));
        var candidates=new List<double>();
        for(int i=0;i<pixels.Length;i++) if(pixels[i]) {
            var queue=new Queue<int>(); queue.Enqueue(i); pixels[i]=false;
            var columns=new HashSet<int>(); var centerYs=new List<int>();
            while(queue.Count>0) {
                int current=queue.Dequeue(),x=current%width,y=current/width;
                columns.Add(x); if(x==30) centerYs.Add(y+top);
                for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++) {
                    int nx=x+dx,ny=y+dy;
                    if(nx>=0 && nx<width && ny>=0 && ny<height) {
                        int n=ny*width+nx;
                        if(pixels[n]) { pixels[n]=false; queue.Enqueue(n); }
                    }
                }
            }
            if(columns.Count>=55 && centerYs.Count>0) {
                centerYs.Sort(); candidates.Add(centerYs[centerYs.Count/2]);
            }
        }
        return candidates.Count==1 ? candidates[0] : (double?)null;
    }
    static string F(double x) { return x.ToString("G17",CultureInfo.InvariantCulture); }
    public static void Run(string imageDirectory,string outputDirectory) {
        Directory.CreateDirectory(outputDirectory);
        string pointsPath=Path.Combine(outputDirectory,"points.csv"), axesPath=Path.Combine(outputDirectory,"calibration.csv");
        if(File.Exists(pointsPath)||File.Exists(axesPath)) throw new IOException("Refusing to overwrite digitized data");
        string[] trials={"camu","box1","box2"};
        string[] states={"G","F","E","L","A","Xy","Xl","Xa"};
        int[][] lefts={new[]{1960,3552,5143,6735},new[]{1485,3076,4668,6260},new[]{1920,3512,5103,6695}};
        int[][] tops={new[]{7456,8993},new[]{6197,7733},new[]{4856,6392}};
        int[][] bottoms={new[]{8649,10185},new[]{7389,8925},new[]{6049,7585}};
        double[][] maxima={new[]{60.0,60.0,25.0,10.0,6.0,1.0,0.6,0.002},new[]{60.0,60.0,7.0,5.0,16.0,0.3,0.5,0.010},new[]{50.0,70.0,8.0,5.0,20.0,0.3,0.75,0.030}};
        // Stay inside the right axis so the window does not clip at t=144.
        int[] times={12,24,36,48,72,96,120};
        using(var points=new StreamWriter(pointsPath)) using(var axes=new StreamWriter(axesPath)) {
            points.WriteLine("experiment,variable,time_h,value,unit,pixel_x,pixel_y,axis_span,thesis_page,figure");
            axes.WriteLine("experiment,variable,x_at_0h,x_at_144h,y_at_0,y_at_labeled_max,labeled_max,plot_top,plot_bottom,thesis_page");
            for(int trial=0;trial<3;trial++) using(var b=new Bitmap(Path.Combine(imageDirectory,"cocoa-thesis-page-"+(45+trial)+".png"))) {
                if(b.Width!=9000 || b.Height!=12729) throw new ArgumentException("Expected a 9000 x 12729 rendering");
                for(int state=0;state<8;state++) {
                    int left=lefts[trial][state%4],top=tops[trial][state/4],bottom=bottoms[trial][state/4];
                    var darkX=new List<int>(); var darkY=new List<int>();
                    for(int x=left-15;x<left+1100;x++) if(Dark(b.GetPixel(x,top-15))) darkX.Add(x);
                    for(int y=top-10;y<=bottom+10;y++) if(Dark(b.GetPixel(left-15,y))) darkY.Add(y);
                    var xticks=Centers(darkX); var yticks=Centers(darkY);
                    if(xticks.Count!=7 || yticks.Count<3) throw new ArgumentException("Unexpected ticks: "+trials[trial]+" "+states[state]+" x="+xticks.Count+" y="+yticks.Count);
                    double x0=xticks[0],x144=xticks[6],ymax=yticks[0],y0=yticks[yticks.Count-1],maximum=maxima[trial][state];
                    axes.WriteLine(string.Join(",",trials[trial],states[state],F(x0),F(x144),F(y0),F(ymax),F(maximum),top,bottom,45+trial));
                    foreach(int time in times) {
                        int x=(int)Math.Round(x0+(x144-x0)*time/144.0);
                        double? y=SolidY(b,x,top+3,bottom-3);
                        if(!y.HasValue) { Console.WriteLine("Omitted occluded/ambiguous point: "+trials[trial]+" "+states[state]+" "+time); continue; }
                        double value=(y0-y.Value)/(y0-ymax)*maximum;
                        points.WriteLine(string.Join(",",trials[trial],states[state],time,F(value),"mg/g pulp",x,F(y.Value),F(maximum),45+trial,"2."+(3+trial)));
                    }
                }
            }
        }
    }
}
'@
[CocoaFigureDigitizer]::Run((Resolve-Path -LiteralPath $ImageDirectory).Path,$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OutputDirectory))
