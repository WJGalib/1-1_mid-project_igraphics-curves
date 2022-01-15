#include<stdio.h>
#include<stdlib.h>
#include"iGraphics.h"
#include<math.h>
#include<string.h>
#define PI 3.141592653589793238462643383279502884197169
#define MATTE 0x00aaaaaa

int	n=3, pathAnimState=0, h1=0, drag2=0, v1=0, v2=0, cont_mode=0, dot=0, drag=0,
MB=799, colIn[6]={0,0,24,86,1,1}, tracer_img[100]={0}, tweakIn[5]={0,0, 24, iScreenHeight-224, 1}, 
paused[100]={0}, path_hid[100]={0}, tracer_hid[100]={0}, dir[100], 
blinkAnimState=0, highlightSum=0, textIn=0, help=0, delIn[3]={0,0,0}, origin[2]={680,350}, 
origin2[2]={680,350}, f_mode=0, type_G=0, deg_P=1;

// strings of numerical variables are needed to display with iText()

char nStr[3]="03", tweakNoStr[20]="01", ampStr[20], trcVisStr[20]="   Hide",
frStr[20], elevStr[20], phaseStr[20], visStr[20]="      Hide", sincosStr[20]=" Change to cos", 
imgNames[11][15], eqnDispStr[100], ampDispStr[20], frDispStr[20], phaseDispStr[20], 
elevDispStr[20], sincosDispStr[5], G_ampStr[10], G_frStr[10], type_GStr[21][10], aStr[20], 
pStr[20], n_PStr[20], deg_PStr[20], b_GStr[20], c_arrStr[130], fStr[250], dxStr[20]="5.00",
col1Str[4]="000", col2Str[4]="000", col3Str[4]="000";

double a=1, p=1, n_P=1, b_G=2.718, c_arr[13]={2, 0.5};

// Each index [n] of arrays of properties of curves represents the properties for nth curve
// index [0] may represent property o the sum, some other entity using the property, or nothing at all
// for color value arrays, [0],[100],[101],[102] represent the colour of the sum, bg, grid, and axes respectively

double G_amp=1, G_fr=1, amp[100]={0, 60, 180, 80}, fr[100]={0,1,2,3}, dx=5, y_i[100]={0}, y_ip1[100]={0}, 
zoomFac[2]={1,1}, sinosum=0, amp2[100]={0, 60, 180, 80}, fr2[100]={0,1,2,3}, 
tracer_x[100]={50, 150, 320, 0}, unit=100, elev[100]={0}, elev2[100]={0},
tracer_phase[100]={0, 150, 320, 50}, phase[100]={0}, sincos[100]={0, 0, 0, 90}, phase2[100]={0}, 
col_r[103]={220,255,87,136}, col_g[103]={220,98,208,190}, col_b[103]={220,102,108,210}, 
col_h[102]={0}, col_s[102]={0}, col_v[102]={0};

double arrsum (double x[]) {
    int i; double sum=0;
    for (i=1; i<=n; i++) sum+=x[i];
    return sum; 
};

double sgn (double x) {
    if (x>=0) return 1;
    else return -1;
};

// General functions for representing each sin/cos surve

double sinusoid (int var, int j) {
    return  
    (  amp[j] * G_amp*zoomFac[1] * sin( (G_fr*fr[j]* (var- origin[0])/zoomFac[0] + phase[j] + sincos[j]) *PI/180 ) + elev[j]*zoomFac[1] + origin[1] );
};

double d_dxSinusoid (int var, int j) {
    return  
    (  amp[j] * G_amp * G_fr*fr[j]*PI/180 * cos( (G_fr*fr[j]* (var- origin[0])/zoomFac[0] + phase[j] + sincos[j]) *PI/180 )  );
};

double inverseSinusoid (int var, int j) {
    return 
    (   ( 180/PI * asin ((var-origin[1]-elev[j]*zoomFac[1])/(amp[j]*G_amp*zoomFac[1])) - phase[j] - sincos[j] ) * zoomFac[0]/G_fr/fr[j] + origin[0]   );
};

// MATH FOR FOURIER TRANSFORM---

double P (double x, double n, int deg, double c[]) {
    int i; double sum=0;
    for (i=0; i<=deg; i++) sum += ( c[i] * pow(x, i) );
    return (  pow( sum, n) );
};

double G (double x, int type, double b) {
    if (type==0)       return (  x  );
    else if (type==1)  return (  sin (x)  );
    else if (type==2)  return (  cos (x)  );
    else if (type==3)  return (  tan (x)  );
    else if (type==4)  return (  1/sin(x)  );
    else if (type==5)  return (  1/cos(x)  );
    else if (type==6)  return (  1/tan (x)  );
    else if (type==7)  return (  asin (x)  );
    else if (type==8)  return (  acos (x)  );
    else if (type==9)  return (  atan (x)  );
    else if (type==10) return (  asin (1/x)  );
    else if (type==11) return (  acos (1/x)  );
    else if (type==12) return (  PI/2 - atan (x)  );
    else if (type==13) return (  sinh (x)  );
    else if (type==14) return (  cosh (x)  );
    else if (type==15) return (  tanh (x)  );
    else if (type==16) return (  1/sinh (x)  );
    else if (type==17) return (  1/cosh (x)  );
    else if (type==18) return (  1/tanh (x)  );
    else if (type==19) return (  log (x)/log(b)  );
    else if (type==20) return (  pow(b, (x))  );
};

double f (double x, int type_f, int type_G, double n_P, int deg_P, double c_P[], double b_G, double a, double p, double m, double L) {    
    if (type_f == 0)      return  (  cos(2*PI*x*m/L) * a * pow( G ( P(x, n_P, deg_P, c_P), type_G, b_G ), p )  );  
    else if (type_f == 1) return  (  sin(2*PI*x*m/L) * a * pow( G ( P(x, n_P, deg_P, c_P), type_G, b_G ), p )  );  
};

double integral (double L1, double L2, int type_f, int type_G, double n_P, int deg_P, double c_P[], double b_G, double a, double p, double m) {
    double i, value=0;
    for (i=L1; i<=L2; i+=0.001) {
        value += 0.001 * f(i, type_f, type_G, n_P, deg_P, c_P, b_G, a, p, m, (L2-L1));
    };
    return value;
};


// FUNCTIONS FOR ADDING AND DELETING CURVES:

void deleteCurve (int x) {
    int i;
    for (i=x; i<n; i++) {
        path_hid[i]= path_hid[i+1]; tracer_hid[i] = tracer_hid[i+1]; dir[i] = dir[i+1]; paused[i] = paused[i+1];
        amp[i] = amp[i+1]; fr[i] = fr[i+1]; y_i[i] = y_i[i+1]; y_ip1[i] = y_ip1[i+1]; amp2[i] = amp2[i+1];
        fr2[i] = fr2[i+1]; tracer_x[i] = tracer_x[i+1]; col_r[i] = col_r[i+1]; col_g[i] = col_g[i+1]; 
        col_b[i] = col_b[i+1]; phase[i] = phase[i+1]; sincos[i] = sincos[i+1]; phase2[i] = phase2[i+1]; 
        elev[i] = elev[i+1]; elev2[i] = elev2[i+1]; tracer_img[i] = tracer_img[i+1]; tracer_phase[i] = tracer_phase[i+1];
    };
    if (x==n) {
        tweakIn[4]--;
        tweakNoStr[0]='0'+ tweakIn[4]/10;
        tweakNoStr[1]='0'+ tweakIn[4]%10;
    };
    n--;
    nStr[0]='0'+ n/10;
    nStr[1]='0'+ n%10;
};

void addCurve () {
    n++;
    nStr[0]='0'+ n/10;
    nStr[1]='0'+ n%10;
    path_hid[n] = 0; tracer_hid[n] = 0; dir[n] = 1; paused[n] = 0; sincos[n] = 90*abs(rand()%2); elev[n] = 0; tracer_img[n]=0;
    amp[n] = (45+rand()%150)/zoomFac[1]; fr[n] = (0.5+fabs(rand()%20/10.0))*zoomFac[0]; tracer_phase[n] = tracer_x[n] = (27+fabs(rand()%40));
    phase[n] = (30+rand()%55)*zoomFac[0]; col_r[n] = 120+rand()%100; col_g[n] = 100+rand()%90; col_b[n] = 101+rand()%100;
    tweakIn[4] = n;
    tweakIn[0] = 2;
};

// Math for interconversion of colour modes (RGB and HSV)

void hsvrgb (double hsv[], double rgb[]) {
    double H=hsv[0],S=hsv[1],V=hsv[2],R,G,B, C,X,m;
    S/=100; V/=100;
    C = 1.0* V * S;
    H /= 60.0;
    X = 1.0* C * (  1.0 - fabs( fmod(H,2.0) - 1.0 )  );
    m = 1.0*(V-C);
    if (0<=H && H<1) {
        R=C; G=X; B=0;
    } else if (1<=H && H<2) {
        R=X; G=C; B=0;
    } else if (2<=H && H<3) {
        R=0; G=C; B=X;
    } else if (3<=H && H<4) {
        R=0; G=X; B=C;
    } else if (4<=H && H<5) {
        R=X; G=0; B=C;
    } else if (5<=H && H<6) {
        R=C; G=0; B=X;
    };
    rgb[0] = (R+m)*255;
    rgb[1] = (G+m)*255;
    rgb[2] = (B+m)*255;
};

double hsvrgb_R (double H, double S, double V) {
    double R;
    S/=100; V/=100;
    double C = V * S;
    H /= 60.0;
    double X = C * (  1.0 - fabs( fmod(H,2.0) - 1.0 )  );
    double m = 1.0*(V-C);
    if (0<=H && H<1)      R=C;
    else if (1<=H && H<2) R=X; 
    else if (2<=H && H<3) R=0;
    else if (3<=H && H<4) R=0;
    else if (4<=H && H<5) R=X;
    else if (5<=H && H<6) R=C;
    return ( (R+m)*255 );
};

double hsvrgb_G (double H, double S, double V) {
    double G;
    S/=100; V/=100;
    double C = V * S;
    H /= 60.0;
    double X = C * (  1.0 - fabs( fmod(H,2.0) - 1.0 )  );
    double m = 1.0*(V-C);
    if (0<=H && H<1)      G=X;
    else if (1<=H && H<2) G=C; 
    else if (2<=H && H<3) G=C;
    else if (3<=H && H<4) G=X;
    else if (4<=H && H<5) G=0;
    else if (5<=H && H<6) G=0;
    return ( (G+m)*255 );
};

double hsvrgb_B (double H, double S, double V) {
    double B;
    S/=100; V/=100;
    double C = V * S;
    H /= 60.0;
    double X = C * (  1.0 - fabs( fmod(H,2.0) - 1.0 )  );
    double m = 1.0*(V-C);
    if (0<=H && H<1)      B=0;
    else if (1<=H && H<2) B=0; 
    else if (2<=H && H<3) B=X;
    else if (3<=H && H<4) B=C;
    else if (4<=H && H<5) B=C;
    else if (5<=H && H<6) B=X;
    return ( (B+m)*255 );
};

void rgbhsv (double rgb[], double hsv[])  {
    double H,S,V,R=rgb[0],G=rgb[1],B=rgb[2], C;
    R/=255; G/=255; B/=255;
    if (R<=G) {
        if (R<=B) C = R;
        else C = B;
    } else {
        if (G<=B) C = G; 
        else C = B;
    };
    if (R>=G) {
        if (R>=B) V = R;
        else V = B;
    } else {
        if (G>=B) V = G;
        else V = B;
    };
    if (V==C) H = 0;
    else if (V==R) H = 60.0 * ( ((G-B)/(V-C)) + 0 );
    else if (V==G) H = 60.0 * ( ((B-R)/(V-C)) + 2 );
    else if (V==B) H = 60.0 * ( ((R-G)/(V-C)) + 4 );

    if (V==0) S = 0;
    else S = (V-C)/V;
    if (H<0) H+=360;
    if (H>359) H-=359;
    hsv[0] = H;
    hsv[2] = V*100; 
    hsv[1] = S*100;
};

// Function to update the display strings in colour menu
void colDispInit (int p) {
    double rgb1[6], hsv1[6];
    if (colIn[5]) {
        rgb1[0] = col_r[p]; rgb1[1] = col_g[p]; rgb1[2] = col_b[p];
        rgbhsv(rgb1, hsv1);
        col_h[p]=hsv1[0]; col_s[p]=hsv1[1]; col_v[p]=hsv1[2];
        col3Str[0] = '0' + (int)col_v[p]/100;
        col3Str[1] = '0' + ((int)col_v[p]/10)%10;
        col3Str[2] = '0' + (int)col_v[p]%10;
        col2Str[0] = '0' + (int)col_s[p]/100;
        col2Str[1] = '0' + ((int)col_s[p]/10)%10;
        col2Str[2] = '0' + (int)col_s[p]%10;
        col1Str[0] = '0' + (int)col_h[p]/100;
        col1Str[1] = '0' + ((int)col_h[p]/10)%10;
        col1Str[2] = '0' + (int)col_h[p]%10;
    } else {
        col3Str[0] = '0' + (int)col_b[p]/100;
        col3Str[1] = '0' + ((int)col_b[p]/10)%10;
        col3Str[2] = '0' + (int)col_b[p]%10;
        col2Str[0] = '0' + (int)col_g[p]/100;
        col2Str[1] = '0' + ((int)col_g[p]/10)%10;
        col2Str[2] = '0' + (int)col_g[p]%10;
        col1Str[0] = '0' + (int)col_r[p]/100;
        col1Str[1] = '0' + ((int)col_r[p]/10)%10;
        col1Str[2] = '0' + (int)col_r[p]%10;
    };
};

// Function to update display strings in edit/add menu
void tweakDispInit() {
    nStr[0]='0'+ n/10;
    nStr[1]='0'+ n%10;
    snprintf (phaseStr, 10, "%lg", phase[tweakIn[4]]);
    snprintf (ampStr, 10, "%lg", amp[tweakIn[4]]);
    snprintf (frStr, 10, "%lg", fr[tweakIn[4]]);
    snprintf (elevStr, 10, "%lg", elev[tweakIn[4]]);
    snprintf (tweakNoStr, 10, "%02d", tweakIn[4]);
    if ( path_hid[tweakIn[4]] ) strcpy(visStr, "     Show");
    else strcpy(visStr, "      Hide");
    if ( sincos[tweakIn[4]] ) strcpy(sincosStr, " Change to sin");
    else strcpy(sincosStr, " Change to cos");
    if ( tracer_hid[tweakIn[4]] ) strcpy(trcVisStr, "  Show");
    else strcpy(trcVisStr, "   Hide");
};

//Function to automatically change sum/grid color based on bgcolor and sum display mode
void updateColTheme() {
    double rgb1[6], hsv1[6];
    if (!colIn[5]) {
        rgb1[0] = col_r[101]; rgb1[1] = col_g[101]; rgb1[2] = col_b[101];
        rgbhsv(rgb1, hsv1);
        col_h[101]=hsv1[0]; col_s[101]=hsv1[1]; col_v[101]=hsv1[2];
    };
   
    if (col_v[101]>=50) {
        col_r[100] = col_r[101] * 0.85;
        col_g[100] = col_g[101] * 0.85;
        col_b[100] = col_b[101] * 0.85;
        if (highlightSum) {
            col_r[0]=50; col_g[0]=50; col_b[0]=50;
            col_r[102]=col_r[100]; col_g[102]=col_g[100]; col_b[102]=col_b[100];
        } else {
           col_r[0]=col_r[100]; col_g[0]=col_g[100]; col_b[0]=col_b[100];
           col_r[102]=50; col_g[102]=50; col_b[102]=50;
        };
    } else {
        col_r[100] = col_r[101] + (255 - col_r[101]) * 0.15;
        col_g[100] = col_g[101] + (255 - col_g[101]) * 0.15;
        col_b[100] = col_b[101] + (255 - col_b[101]) * 0.15;
        if (highlightSum) {
            col_r[0]=200; col_g[0]=200; col_b[0]=200;
            col_r[102]=col_r[100]; col_g[102]=col_g[100]; col_b[102]=col_b[100];
        } else {
            col_r[0]=col_r[100]; col_g[0]=col_g[100]; col_b[0]=col_b[100];
            col_r[102]=200; col_g[102]=200; col_b[102]=200;
        };
    };
};

//Function to display UI buttons all the time
void drawButtons() {
    if (dot) iShowBMP2 (iScreenWidth-72, iScreenHeight-60, "dot.bmp", MATTE);
    else iShowBMP2 (iScreenWidth-72, iScreenHeight-60, "line.bmp", MATTE);
    if (cont_mode) iShowBMP2 (iScreenWidth-132, iScreenHeight-60, "reenter.bmp", MATTE);
    else iShowBMP2 (iScreenWidth-132, iScreenHeight-60, "bounce.bmp", MATTE);
    iShowBMP2 (iScreenWidth-192, iScreenHeight-60, "bgcol.bmp", MATTE);
    iShowBMP2 (24, iScreenHeight-60, "add.bmp", MATTE);
    iShowBMP2 (84, iScreenHeight-60, "edit.bmp", MATTE);
    iShowBMP2 (144, iScreenHeight-60, "tracer.bmp", MATTE);
    iShowBMP2 (204, iScreenHeight-60, "sum.bmp", MATTE);
    iShowBMP2 (264, iScreenHeight-60, "fourier.bmp", MATTE);
    iShowBMP2 (324, iScreenHeight-60, "delete.bmp", MATTE);
    iShowBMP2 (24, 24, "origin.bmp", MATTE);
    iShowBMP2 (84, 24, "helpb.bmp", MATTE);
};

//Function to show delete confirmation popup
void drawDelMenu (int x, int y) {
    iSetColor (40, 49, 71);
    iFilledRectangle (x, y-100, 200, 95);
    iSetColor (220,220,220);
    iText(x+24, y-36, "Delete all curves? Confirm?", GLUT_BITMAP_HELVETICA_12);
    iSetColor (42, 76, 130);
    iRectangle (x, y-100, 200, 95);
    iSetColor (112, 146, 190);
    iFilledRectangle (x+10, y-90, 85, 36);
    iFilledRectangle (x+105, y-90, 85, 36);
    iSetColor (220,220,220);
    iRectangle (x+10, y-90, 85, 36);
    iRectangle (x+105, y-90, 85, 36);
    iSetColor (0,0,0);
    iText(x+40, y-78, "Yes!", GLUT_BITMAP_HELVETICA_12);
    iText(x+125, y-78, "Cancel", GLUT_BITMAP_HELVETICA_12);
};

// Function to display equation of curve when add/edit menu is open, and update relavent strings
void drawEqnDisp() {

    iShowBMP2 (iScreenWidth/2-378, 48, "eqnbgL.bmp",  MATTE);
    iShowBMP2 (iScreenWidth/2-252, 48, "eqnbg.bmp",  MATTE);
    iShowBMP2 (iScreenWidth/2-126, 48, "eqnbg.bmp",  MATTE);
    iShowBMP2 (iScreenWidth/2, 48, "eqnbg.bmp",  MATTE);
    iShowBMP2 (iScreenWidth/2+126, 48, "eqnbg.bmp",  MATTE);
    iShowBMP2 (iScreenWidth/2+252, 48, "eqnbgR.bmp",  MATTE);
    
    if ( sincos[tweakIn[4]] ) strcpy(sincosDispStr, "cos");
    else strncpy(sincosDispStr, "sin", 4);

    if ( amp[tweakIn[4]]==1 ) strcpy(ampDispStr, "");
    else if (amp[tweakIn[4]]==-1) strcpy(ampDispStr, "-");
    else snprintf (ampDispStr, 15, "%lg", amp[tweakIn[4]]);

    if ( phase[tweakIn[4]]==0 ) strcpy(phaseDispStr, "");
    else snprintf (phaseDispStr, 15, "%+lg", phase[tweakIn[4]]*PI/180);

    if (fabs(fr[tweakIn[4]]-180/PI) < 0.001) strcpy(frDispStr, "");
    else if (fabs(fr[tweakIn[4]]+180/PI) < 0.001) strcpy(frDispStr, "-");
    else snprintf (frDispStr, 15, "%lg", fr[tweakIn[4]]*PI/180);

    if ( elev[tweakIn[4]]==0 ) strcpy(elevDispStr, "");
    else snprintf (elevDispStr, 15, "%+lg", elev[tweakIn[4]]);

    iSetColor (220,220,220);
    if (amp[tweakIn[4]] == 0) {
        if (elev[tweakIn[4]]==0) snprintf(eqnDispStr, 99, "y = 0", elevDispStr);
        else snprintf(eqnDispStr, 99, "y = %s", elevDispStr);
    } else if (fr[tweakIn[4]] == 0) snprintf (eqnDispStr, 99, "y = %s %s (%s) %s", ampDispStr, sincosDispStr, phaseDispStr, elevDispStr);
    else snprintf (eqnDispStr, 99, "y = %s %s (%sx %s) %s", ampDispStr, sincosDispStr, frDispStr, phaseDispStr, elevDispStr);
    iText (iScreenWidth/2-10*strlen(eqnDispStr)/2, 78, eqnDispStr, GLUT_BITMAP_TIMES_ROMAN_24);
};

// Function to draw pacman if it is set as tracer 
void drawPacman (int state, int x, int y, int j) {
    iSetColor(249,237,0);
    double pacmanX[8]={x, x+12, x+16, x+12, x, x-12, x-16, x-12}, 
    pacmanY[8]={y+16, y+12, y, y-12, y-16, y-12, y, y+12};
    if (dir[j] > 0) {
        if (state) {
            pacmanX[2] = x;
            pacmanY[2] = y-8;
        };
        iFilledPolygon (pacmanX, pacmanY, 8);
        iSetColor(0,0,0);
        if (!state) iLine(x+4,y,x+16,y);
    } else {
        if (state) {
            pacmanX[6] = x;
            pacmanY[6] = y-8;
        };
        iFilledPolygon (pacmanX, pacmanY, 8);
        iSetColor(0,0,0);
        if (!state) iLine(x-4,y,x-16,y);
    }
    iPoint (x, y+8, 2);
};

// Function to display global amp/freq factors
void drawGFacs() {
    iShowBMP2 (iScreenWidth-200, 0, "gfacbg.bmp", MATTE);
    iSetColor (220,220,220);
    iText (iScreenWidth-178, 40, "Global amplitude factor:", GLUT_BITMAP_HELVETICA_10);
    if (G_amp>0) iText (iScreenWidth-60, 40, G_ampStr, GLUT_BITMAP_HELVETICA_10);
    else iText (iScreenWidth-60, 40, "0", GLUT_BITMAP_HELVETICA_10);    
    iText (iScreenWidth-178, 25, "Global frequency factor:", GLUT_BITMAP_HELVETICA_10);
    if (G_fr>0) iText (iScreenWidth-60, 25, G_frStr, GLUT_BITMAP_HELVETICA_10);
    else iText (iScreenWidth-60, 25, "0", GLUT_BITMAP_HELVETICA_10);
    iText (iScreenWidth-178, 10, "Tracer speed:", GLUT_BITMAP_HELVETICA_10);
    iText (iScreenWidth-88, 10, "-", GLUT_BITMAP_HELVETICA_10);
    iText (iScreenWidth-74, 10, dxStr, GLUT_BITMAP_HELVETICA_10);
    iText (iScreenWidth-42, 10, "+", GLUT_BITMAP_HELVETICA_10);
};

// Function to update strings in the fourier series input menu
void fourierDispInit() {
    int i, j=0; char temp[30], P_str[150]="", aDispStr[20], type_GDispStr[20];
   
    snprintf (aStr, 18, "%lg", a);
    snprintf (pStr, 18, "%lg", p);
    snprintf (n_PStr, 18, "%lg", n_P);
    snprintf (b_GStr, 18, "%lg", b_G);
    snprintf (deg_PStr, 18, "%d", deg_P);
    strncpy (c_arrStr, "", 128);
    for(i=0; i<=deg_P; i++) {
        snprintf (temp, 18, "%lg ", c_arr[i]);
        strncat (c_arrStr, temp, 128);
        if (fabs(c_arr[i])>0.000001) {
            j++;
            if (i==0) {
                if (fabs(fabs(c_arr[i])-1)<0.00001) snprintf (temp, 28, "%s 1 ", ((sgn(c_arr[i])>0)?((j==1)?"":"+"):"-") );
                else snprintf (temp, 28, "%s %lg ", ((sgn(c_arr[i])>0)?((j==1)?"":"+"):"-") , fabs(c_arr[i]) );
            } else if (i==1) {
                if (fabs(fabs(c_arr[i])-1)<0.00001) snprintf (temp, 28, "%s x ", ((sgn(c_arr[i])>0)?((j==1)?"":"+"):"-") );
                else snprintf (temp, 28, "%s %lgx ", ((sgn(c_arr[i])>0)?((j==1)?"":"+"):"-") , fabs(c_arr[i]) );
            } else {
                if (fabs(fabs(c_arr[i])-1)<0.00001) snprintf (temp, 28, "%s x^%d ", ((sgn(c_arr[i])>0)?((j==1)?"":"+"):"-"), i );
                else snprintf (temp, 28, "%s %lgx^%d ", ((sgn(c_arr[i])>0)?((j==1)?"":"+"):"-") , fabs(c_arr[i]), i );
            };
            strncat (P_str, temp, 148);
        };
    };

    if (fabs(a-1)<0.00001) strcpy(aDispStr, "");
    else strcpy(aDispStr, aStr);

    if (type_G==0) strcpy(type_GDispStr, "");
    else if (type_G==19) snprintf(type_GDispStr, 18, "log_(%s)", b_GStr);
    else if (type_G==20) snprintf(type_GDispStr, 18, "(%s)^", b_GStr);
    else strcpy(type_GDispStr, type_GStr[type_G]);

    if (fabs(p-1)<0.00001) {
        if (fabs(n_P-1)<0.00001) snprintf (fStr, 248, "%s %s( %s)", aDispStr, type_GDispStr, P_str);
        else snprintf (fStr, 248, "%s (%s(( %s)^%s))", aDispStr, type_GDispStr, P_str, n_PStr);
    } else if (fabs(n_P-1)<0.00001) {
        snprintf (fStr, 248, "%s (%s( %s))^%s", aDispStr, type_GDispStr, P_str, pStr);
    } else snprintf (fStr, 248, "%s (%s(( %s)^%s))^%s", aDispStr, type_GDispStr, P_str, n_PStr, pStr);

};

//Function to display Fourier series input menu
void drawFourierHelp() {

    iSetColor (28,33,47);
    iFilledRectangle (0,0, 1360, 700);
    iShowBMP(iScreenWidth/2-486, iScreenHeight/2-208, "fourier_help.bmp");
    
    iSetColor (220,220,220);
    iFilledRectangle (220, 220, 74,24);
    if (textIn==6 && blinkAnimState) iSetColor (249,237,101);
    iFilledRectangle (340, 220, 54,24);
    if (textIn==7 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle (448, 220, 54,24);
    if (textIn==8 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle (588, 220, 54,24);
    if (textIn==9 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle (850, 220, 54,24);
    if (textIn==10 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle (788, 178, 380,24);
    if (textIn==11 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle (1018, 220, 54,24);

    double x[3]={224, 234, 234}, y[3]={231, 235, 226};
    iSetColor (0,0,0);
    iFilledPolygon (x,y,3);
    x[0]=290; x[1]=x[2]=280; y[0]=231; y[1]=236; y[2]=226;
    iFilledPolygon (x,y,3);

    iText (240, 226, type_GStr[type_G], GLUT_BITMAP_HELVETICA_10);
    iText (346, 226, aStr, GLUT_BITMAP_HELVETICA_10);
    iText (454, 226, pStr, GLUT_BITMAP_HELVETICA_10);
    iText (594, 226, n_PStr, GLUT_BITMAP_HELVETICA_10);
    iText (856, 226, deg_PStr, GLUT_BITMAP_HELVETICA_10);
    iText (794, 184, c_arrStr, GLUT_BITMAP_HELVETICA_10);
    iText (1024, 226, b_GStr, GLUT_BITMAP_HELVETICA_10);
    iSetColor (220,220,220);
    iText (iScreenWidth/2-8.5*strlen(fStr)/2, 126, fStr, GLUT_BITMAP_9_BY_15);
    iSetColor (40, 49, 71);
    iFilledCircle (iScreenWidth-40, iScreenHeight-40, 16);
    iFilledRectangle (iScreenWidth/2-32, 24, 76, 36);
    iSetColor (220,220,220);
    iCircle (iScreenWidth-40, iScreenHeight-40, 16);
    iRectangle (iScreenWidth/2-32, 24, 76, 36);
    iText (iScreenWidth/2-24, 35, "PLOT!", GLUT_BITMAP_HELVETICA_18);
    iText (iScreenWidth-47, iScreenHeight-46, "X", GLUT_BITMAP_HELVETICA_18);
    
};


void drawHelp() {
    iSetColor (28,33,47);
    iFilledRectangle (0,0, 1360, 700);
};

// Function to display add/edit menu
void drawTweakMenu (int x, int y) {
    iSetColor (40, 49, 71);
    iFilledRectangle ( x, y, 340, 200);
    double arrowX[3] = {x+68, x+82,x+82}, arrowY[3] = {y-25, y-30, y-20};

    iFilledRectangle ( x, y-72, 340, 66);
    iSetColor (col_r[tweakIn[4]], col_g[tweakIn[4]], col_b[tweakIn[4]]);
    if (tracer_img[tweakIn[4]]==0) iFilledCircle ( x+307, y-39, 6 );
    else if (tracer_img[tweakIn[4]]==1) iFilledRectangle ( x+299, y-47, 16, 16);
    else if (tracer_img[tweakIn[4]]==10) drawPacman (blinkAnimState, x+307, y-39, 100);
    else if (tracer_img[tweakIn[4]]==9) iShowBMP2 (x+282, y-64, "d_dx.bmp", MATTE);
    else iShowBMP2 (x+282, y-64, imgNames[tracer_img[tweakIn[4]]], MATTE);
    
    iSetColor (220, 220, 220);
    iText (x+8, y-30, "Tracer: ", GLUT_BITMAP_HELVETICA_12);
    iFilledPolygon (arrowX, arrowY, 3);
    iText (x+88, y-30, imgNames[tracer_img[tweakIn[4]]], GLUT_BITMAP_HELVETICA_12);
    arrowX[0]=x+187; arrowX[1]=x+173; arrowX[2]=x+173;
    iFilledPolygon (arrowX, arrowY, 3);
    iText (x+8, y-61, "Leading", GLUT_BITMAP_HELVETICA_12);
    iFilledRectangle (x+82, y-57, 186, 3);
    
    iSetColor (150,150,150);
    iFilledRectangle ( x+80 + tracer_phase[tweakIn[4]]/7.3118, y-67, 5, 24 );

    iSetColor (42, 76, 130);
    iRectangle ( x, y, 340, 200);
    iRectangle ( x, y-72, 340, 66);

    iSetColor (112, 146, 190);
    iFilledRectangle ( x+218, y-38, 50, 24 );

    iFilledRectangle ( x+257, y+8, 75, 24 );
    iFilledRectangle ( x+174, y+8, 75, 24 );
    iFilledRectangle ( x+8, y+8, 95, 24 );
    if (tweakIn[0]==2) iSetColor (102, 225, 139);
    else iSetColor (250,107,107);
    iFilledRectangle ( x+111, y+8, 55, 24 );
    

    iSetColor (220, 220, 220);
    iRectangle ( x+218, y-38, 50, 24 );

    iRectangle ( x+257, y+8, 75, 24 );
    iRectangle ( x+174, y+8, 75, 24 );
    iRectangle ( x+111, y+8, 55, 24 );
    iRectangle ( x+8, y+8, 95, 24 );

    iSetColor (0,0,0);
    iText (x+218, y-30, trcVisStr, GLUT_BITMAP_HELVETICA_12);
    iText (x+10, y+15, sincosStr, GLUT_BITMAP_HELVETICA_12);
    if (tweakIn[0]!=2) iText (x+91, y+15, "       Delete", GLUT_BITMAP_HELVETICA_12);
    else iText (x+91, y+15, "        Add", GLUT_BITMAP_HELVETICA_12);
    iText (x+174, y+15, visStr, GLUT_BITMAP_HELVETICA_12);
    iText (x+257, y+15, "     Color", GLUT_BITMAP_HELVETICA_12);

    iSetColor (220, 220, 220);
    iText (x+8, y+46, "Phase", GLUT_BITMAP_HELVETICA_12);
    iText (x+8, y+78, "Elevation", GLUT_BITMAP_HELVETICA_12);
    iText (x+8, y+110, "Frequency", GLUT_BITMAP_HELVETICA_12);
    iText (x+8, y+142, "Amplitude", GLUT_BITMAP_HELVETICA_12);
    iText (x+8, y+174, "Curve no.", GLUT_BITMAP_HELVETICA_12);

    if (tweakIn[0]!=2) {
        arrowY[0]=y+179; arrowY[1]=y+184; arrowY[2]=y+174;
        arrowX[0]=x+110; arrowX[1]=arrowX[2]=x+124;
        iFilledPolygon (arrowX, arrowY, 3);
        iText (x+140, y+174, tweakNoStr, GLUT_BITMAP_HELVETICA_12);
        iText (x+170, y+174, "of ", GLUT_BITMAP_HELVETICA_12);
        iText (x+200, y+174, nStr, GLUT_BITMAP_HELVETICA_12);
        arrowX[0]=x+249; arrowX[1]=arrowX[2]=x+235;
        iFilledPolygon (arrowX, arrowY, 3);
    } else {
        iText (x+170, y+174, nStr, GLUT_BITMAP_HELVETICA_12);
        iText (x+185, y+174, "(New)", GLUT_BITMAP_HELVETICA_12);
    };

    if (textIn==5 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle ( x+278, y+40, 54, 24 );
    if (textIn==4 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle ( x+278, y+72, 54, 24 );
    if (textIn==3 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle ( x+278, y+104, 54, 24 );
    if (textIn==2 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    iFilledRectangle ( x+278, y+136, 54, 24 );
    if (textIn==1 && blinkAnimState) iSetColor (249,237,101);
    else iSetColor (220,220,220);
    if (tweakIn[0]!=2) iFilledRectangle ( x+278, y+168, 54, 24 );
    iSetColor (0,0,0);
    iRectangle ( x+278, y+40, 54, 24 );
    iRectangle ( x+278, y+72, 54, 24 );
    iRectangle ( x+278, y+104, 54, 24 );
    iRectangle ( x+278, y+136, 54, 24 );
    if (tweakIn[0]!=2) iRectangle ( x+278, y+168, 54, 24 );

    iText (x+282, y+46, phaseStr, GLUT_BITMAP_HELVETICA_10);
    iText (x+282, y+78, elevStr, GLUT_BITMAP_HELVETICA_10);
    iText (x+282, y+110, frStr, GLUT_BITMAP_HELVETICA_10);
    iText (x+282, y+142, ampStr, GLUT_BITMAP_HELVETICA_10);
    if (tweakIn[0]!=2) iText (x+284, y+174, tweakNoStr, GLUT_BITMAP_HELVETICA_10);

    iSetColor (220,220,220);
    iFilledRectangle (x+82, y+50, 186, 3);
    iFilledRectangle (x+82, y+82, 186, 3);
    iFilledRectangle (x+82, y+114, 186, 3);
    iFilledRectangle (x+82, y+146, 186, 3);

    iSetColor (150,150,150);
    iFilledRectangle ( x+80 + 93 + (fmod(phase[tweakIn[4]],360.0))/3.871, y+40, 5, 24 );
    iFilledRectangle ( x+80 + 93 + ( (fabs(elev[tweakIn[4]])<473)? elev[tweakIn[4]]/5.086 : 93*sgn(elev[tweakIn[4]]) ), y+72, 5, 24 );
    iFilledRectangle ( x+80 + 93 + ( (fabs(fr[tweakIn[4]])<6)? fr[tweakIn[4]]*15.5 : 93*sgn(fr[tweakIn[4]])), y+104, 5, 24 );
    iFilledRectangle ( x+80 + 93 + ( (fabs(amp[tweakIn[4]]) < 472.95)? amp[tweakIn[4]]/5.086 : 93*sgn(amp[tweakIn[4]]) ), y+136, 5, 24 );

    iSetColor (0,0,0);
    iRectangle ( x+80 + 93 + (fmod(phase[tweakIn[4]],360.0))/3.871, y+40, 5, 24 );
    iRectangle ( x+80 + 93 + ( (fabs(elev[tweakIn[4]])<473)? elev[tweakIn[4]]/5.086 : 93*sgn(elev[tweakIn[4]]) ), y+72, 5, 24 );
    iRectangle ( x+80 + 93 + ( (fabs(fr[tweakIn[4]])<6)? fr[tweakIn[4]]*15.5 : 93*sgn(fr[tweakIn[4]])), y+104, 5, 24 );
    iRectangle ( x+80 + 93 + ( (fabs(amp[tweakIn[4]]) < 472.95)? amp[tweakIn[4]]/5.086 : 93*sgn(amp[tweakIn[4]]) ), y+136, 5, 24 );
    
    iRectangle ( x+80 + tracer_phase[tweakIn[4]]/7.3118, y-67, 5, 24 );

    drawEqnDisp();

};

// Function to draw a colour input menu
void drawColMenu (int x, int y) {

    int i, j, k;
    double rgb1[3]={0}, hsv1[3]={0}, hsv2[3]={0};
    iSetColor (40, 49, 71);
    iFilledRectangle (x, y, 340, 160);

    iSetColor (col_r[colIn[4]], col_g[colIn[4]], col_b[colIn[4]]);
    iFilledRectangle (x+51, y+8, 281, 48);
    iSetColor (0,0,0);
    iRectangle (x+51, y+8, 281, 48);
    iSetColor (42, 76, 130);
    iRectangle (x, y, 340, 160);
    if (colIn[5]) iFilledRectangle (x+8, y+8, 35, 20);
    else iFilledRectangle (x+8, y+36, 35, 20);
    iSetColor (0,0,0);
    iRectangle (x+8, y+8, 35, 20);
    iRectangle (x+8, y+36, 35, 20);
    iSetColor (220,220,220);
    iText (x+8, y+14, " HSB",  GLUT_BITMAP_HELVETICA_12);
    iText (x+8, y+42, " RGB", GLUT_BITMAP_HELVETICA_12);
    
    iSetColor (220,220,220);
    iFilledRectangle ( 288+x, 64+y, 44, 24 );
    iFilledRectangle ( 288+x, 96+y, 44, 24 );
    iFilledRectangle ( 288+x, 128+y, 44, 24 );
    iSetColor (0,0,0);
    iRectangle ( 288+x, 64+y, 44, 24 );
    iRectangle ( 288+x, 96+y, 44, 24 );
    iRectangle ( 288+x, 128+y, 44, 24 );

    if (colIn[5]) {
        for (i=0; i<256; i++) {
            
            iSetColor ( hsvrgb_R(1.41137*i, 100, 100), hsvrgb_G(1.41137*i, 100, 100), hsvrgb_B(1.41137*i, 100, 100));
            iLine (25+i+x, 135+y, 25+i+x, 146+y);
            rgb1[0] = col_r[colIn[4]]; rgb1[1] = col_g[colIn[4]]; rgb1[2] = col_b[colIn[4]];
            rgbhsv(rgb1, hsv1);
            col_h[colIn[4]]=hsv1[0]; col_s[colIn[4]]=hsv1[1]; col_v[colIn[4]]=hsv1[2];
           
            iSetColor( hsvrgb_R(col_h[colIn[4]], 0.392157*i, 50), hsvrgb_G(col_h[colIn[4]], 0.392157*i, 50), hsvrgb_B(col_h[colIn[4]], 0.392157*i, 50));
            iLine (25+i+x, 103+y, 25+i+x, 114+y);
          
            iSetColor( hsvrgb_R(col_h[colIn[4]], col_s[colIn[4]], 0.392157*i), hsvrgb_G(col_h[colIn[4]], col_s[colIn[4]], 0.392157*i), hsvrgb_B(col_h[colIn[4]], col_s[colIn[4]], 0.392157*i));
            iLine (25+i+x, 71+y, 25+i+x, 82+y);            
        };

        iSetColor (0,0,0);
        iRectangle ( 24+x, 70+y, 256, 12 );
        iRectangle ( 24+x, 102+y, 256, 12 );
        iRectangle ( 24+x, 134+y, 256, 12 );
        
        iSetColor (150,150,150);
        iFilledRectangle ( 22+col_v[colIn[4]]*2.55+x, 64+y, 5, 24 );
        iFilledRectangle ( 22+col_s[colIn[4]]*2.55+x, 96+y, 5, 24 );
        iFilledRectangle ( 22+col_h[colIn[4]]*0.708333+x, 128+y, 5, 24 );
        iSetColor (0,0,0);
        iRectangle ( 22+col_v[colIn[4]]*2.55+x, 64+y, 5, 24 );
        iRectangle ( 22+col_s[colIn[4]]*2.55+x, 96+y, 5, 24 );
        iRectangle ( 22+col_h[colIn[4]]*0.708333+x, 128+y, 5, 24 );

        iSetColor (220,220,220);
        iText (10+x, 71+y, "B", GLUT_BITMAP_HELVETICA_12);
        iText (10+x, 105+y, "S", GLUT_BITMAP_HELVETICA_12);
        iText (10+x, 137+y, "H", GLUT_BITMAP_HELVETICA_12);
        iSetColor (0,0,0);
        iText (296+x, 71+y, col3Str, GLUT_BITMAP_HELVETICA_12);
        iText (296+x, 105+y, col2Str, GLUT_BITMAP_HELVETICA_12);
        iText (296+x, 137+y, col1Str, GLUT_BITMAP_HELVETICA_12);

    } else {

        for (i=0; i<256; i++) {
            iSetColor (0,0,i);
            iLine (25+i+x, 71+y, 25+i+x, 82+y);
            iSetColor (0,i,0);
            iLine (25+i+x, 103+y, 25+i+x, 114+y);
            iSetColor (i,0,0);
            iLine (25+i+x, 135+y, 25+i+x, 146+y);
        };
        iSetColor (0,0,0);
        iRectangle ( 24+x, 70+y, 256, 12 );
        iRectangle ( 24+x, 102+y, 256, 12 );
        iRectangle ( 24+x, 134+y, 256, 12 );    
        iSetColor (150,150,150);
        iFilledRectangle ( 22+col_b[colIn[4]]+x, 64+y, 5, 24 );
        iFilledRectangle ( 22+col_g[colIn[4]]+x, 96+y, 5, 24 );
        iFilledRectangle ( 22+col_r[colIn[4]]+x, 128+y, 5, 24 );
        iSetColor (0,0,0);
        iRectangle ( 22+col_b[colIn[4]]+x, 64+y, 5, 24 );
        iRectangle ( 22+col_g[colIn[4]]+x, 96+y, 5, 24 );
        iRectangle ( 22+col_r[colIn[4]]+x, 128+y, 5, 24 );

        iSetColor (220,220,220);
        iText (10+x, 71+y, "B", GLUT_BITMAP_HELVETICA_12);
        iText (10+x, 105+y, "G", GLUT_BITMAP_HELVETICA_12);
        iText (10+x, 137+y, "R", GLUT_BITMAP_HELVETICA_12);
        iSetColor (0,0,0);
        iText (296+x, 71+y, col3Str, GLUT_BITMAP_HELVETICA_12);
        iText (296+x, 105+y, col2Str, GLUT_BITMAP_HELVETICA_12);
        iText (296+x, 137+y, col1Str, GLUT_BITMAP_HELVETICA_12);

    };
};

// To rerspond to keyboard/mouse interactions that control the add/edit menu, is called from iMouse() and iKeyboard() 
void tweakValue (int mx, int my, int x, int y, int state) {
    if (!drag2 && state == GLUT_UP) {
        if (mx>=tweakIn[2]+110 && mx<=tweakIn[2]+124 && my>=tweakIn[3]+174 && my<=tweakIn[3]+184) {
            if (tweakIn[4]>1) tweakIn[4]--;
            else if (tweakIn[4]==1) tweakIn[4] = n;            
            if (colIn[0]==2) {
                colIn[4] = tweakIn[4];
                colDispInit(tweakIn[4]);
            };
            tweakDispInit();
        } else if (mx>=tweakIn[2]+235 && mx<=tweakIn[2]+249 && my>=tweakIn[3]+174 && my<=tweakIn[3]+184) {
            if (tweakIn[4]<n) tweakIn[4]++;
            else if (tweakIn[4]==n) tweakIn[4] = 1;            
            if (colIn[0]==2) {
                colIn[4] = tweakIn[4];
                colDispInit(tweakIn[4]);
            };
            tweakDispInit();
        } else if (mx>=tweakIn[2]+174 && mx<=tweakIn[2]+249 && my>=tweakIn[3]+8 && my<=tweakIn[3]+32) {
            if ( path_hid[tweakIn[4]] ) {
                path_hid[tweakIn[4]] = 0;
                strcpy(visStr, "      Hide");
            } else {
                path_hid[tweakIn[4]] = 1;
                strcpy(visStr, "     Show");
            };
        } else if (mx>=tweakIn[2]+8 && mx<=tweakIn[2]+103 && my>=tweakIn[3]+8 && my<=tweakIn[3]+32) {
            if (f_mode==2) f_mode = 0;
            if ( sincos[tweakIn[4]] ) {
                sincos[tweakIn[4]] = 0;
                strcpy(sincosStr, " Change to cos");
            } else {
                sincos[tweakIn[4]] = 90; 
                strcpy(sincosStr, " Change to sin");
            };
        } else if (mx>=tweakIn[2]+111 && mx<=tweakIn[2]+186 && my>=tweakIn[3]+8 && my<=tweakIn[3]+32) {
            if (f_mode==2) f_mode = 0;
            colIn[0] = 0;
            if (tweakIn[0]==1) {
                if (n>1) deleteCurve(tweakIn[4]);
                tweakDispInit();
            } else {
                tweakIn[0] = 0;
            };
        } else if (mx>=tweakIn[2]+257 && mx<=tweakIn[2]+332 && my>=tweakIn[3]+8 && my<=tweakIn[3]+32) {
            if (colIn[0]) {
                if (colIn[0]==1) colIn[0]=0;
                else {                
                    colIn[0] = 0;
                };
            } else {
                colIn[4] = tweakIn[4];
                colDispInit(tweakIn[4]);
                colIn[2] = tweakIn[2]+346; colIn[3] = iScreenHeight-248;
                colIn[0] = 2; 
            };
        } else if (mx>=tweakIn[2]+68 && mx<=tweakIn[2]+82 && my>=tweakIn[3]-30 && my<=tweakIn[3]-20) {
            if ( tracer_img[tweakIn[4]]==0 ) tracer_img[tweakIn[4]] = 10;
            else tracer_img[tweakIn[4]]--;
        } else if (mx>=tweakIn[2]+173 && mx<=tweakIn[2]+187 && my>=tweakIn[3]-30 && my<=tweakIn[3]-20) {
            if ( tracer_img[tweakIn[4]]==10 ) tracer_img[tweakIn[4]] = 0;
            else tracer_img[tweakIn[4]]++;
        } else if (mx>=tweakIn[2]+218 && mx<=tweakIn[2]+268 && my>=tweakIn[3]-38 && my<=tweakIn[3]-14) {
            if ( tracer_hid[tweakIn[4]] ) {
                tracer_hid[tweakIn[4]] = 0;
                strcpy(trcVisStr, "   Hide");
            } else {
                tracer_hid[tweakIn[4]] = 1;
                strcpy(trcVisStr, "  Show");
            };
        };
        if (mx>=tweakIn[2]+278 && mx<=tweakIn[2]+332 && my>=tweakIn[3]+40 && my<=tweakIn[3]+64) {
            textIn = 5;
        } else if (mx>=tweakIn[2]+278 && mx<=tweakIn[2]+332 && my>=tweakIn[3]+72 && my<=tweakIn[3]+96) {
            textIn = 4;
        } else if (mx>=tweakIn[2]+278 && mx<=tweakIn[2]+332 && my>=tweakIn[3]+104 && my<=tweakIn[3]+128) {
            textIn = 3;
        } else if (mx>=tweakIn[2]+278 && mx<=tweakIn[2]+332 && my>=tweakIn[3]+136 && my<=tweakIn[3]+160) {
            textIn = 2;
        } else if (mx>=tweakIn[2]+278 && mx<=tweakIn[2]+332 && my>=tweakIn[3]+174 && my<=tweakIn[3]+184) {
            textIn = 1;
        } else textIn = 0;
    } else if (!drag2) {
        if (  mx >= 80+x && mx <= 80+186+x && ( (my >= 40+y && my <= 63+y) || tweakIn[1]==1 )   ){
            phase[tweakIn[4]] = 1.0*(mx-80-x-93)*3.871;
            tweakIn[1] = 1;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
        } else if (  mx >= 80+x && mx <= 80+186+x && ( (my >= 72+y && my <= 95+y) || tweakIn[1]==2 )   ){
            elev[tweakIn[4]] = 1.0*(mx-80-x-93)*5.086;
            tweakIn[1] = 2;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
        } else if (  mx >= 80+x && mx <= 80+186+x && ( (my >= 104+y && my <= 127+y) || tweakIn[1]==3 )   ){
            fr[tweakIn[4]] = 1.0*(mx-80-x-93)/15.55;
            tweakIn[1] = 3;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
        } else if (  mx >= 80+x && mx <= 80+186+x && ( (my >= 136+y && my <= 159+y) || tweakIn[1]==4 )   ){
            amp[tweakIn[4]] = 1.0*(mx-80-x-93)*5.086;
            tweakIn[1] = 4;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
        } else if (  mx >= 80+x && mx <= 80+186+x && ( (my >= y-67 && my <= y-43) || tweakIn[1]==5 )   ){
            tracer_phase[tweakIn[4]] = 1.0*(mx-80-x)*7.3118;
            tweakIn[1] = 5;
            tweakDispInit();
        };
    };

};

// To rerspond to keyboard/mouse interactions that control a colour input menu, is called from iMouse() and iKeyboard() 
void colValue (int mx, int my, int x, int y, int p, int state) {
    if (!drag2 && state == GLUT_UP) {
        if (mx>=colIn[2]+8 && mx<=colIn[2]+43 && my>=colIn[3]+8 && my<=colIn[3]+28) {
            if (!colIn[5]) colIn[5] = 1;
            colDispInit(p);
        } else if (mx>=colIn[2]+8 && mx<=colIn[2]+43 && my>=colIn[3]+36 && my<=colIn[3]+56) {
            if (colIn[5]) colIn[5] = 0;
            colDispInit(p);
        };
    } else if (!drag2) {
        double rgb2[6], hsv3[6];
        if (  mx >= 25+x && mx <= 25+255+x && ( (my >= 71+y && my <= 82+y) || colIn[1]==1 )   ){
            if (colIn[5]) {
                col_v[p] = (mx-25-x) /2.55;
                hsv3[0] = col_h[p]; hsv3[1] = col_s[p]; hsv3[2] = col_v[p];
                hsvrgb(hsv3, rgb2);
                col_r[p] = rgb2[0]; col_g[p] = rgb2[1]; col_b[p] = rgb2[2];
                col3Str[0] = '0' + (int)col_v[p]/100;
                col3Str[1] = '0' + ((int)col_v[p]/10)%10;
                col3Str[2] = '0' + (int)col_v[p]%10; 
            } else {
                col_b[p] = mx-25-x;
                col3Str[0] = '0' + (int)col_b[p]/100;
                col3Str[1] = '0' + ((int)col_b[p]/10)%10;
                col3Str[2] = '0' + (int)col_b[p]%10; 
            };
            if (colIn[0]==1) updateColTheme();
            colIn[1] = 1;        
        } else if (mx >= 25+x && mx <= 25+255+x && ( (my >= 103+y && my <= 114+y) || colIn[1]==2)  ) {
            if (colIn[5]) {
                col_s[p] = (mx-25-x) /2.55;
                hsv3[0] = col_h[p]; hsv3[1] = col_s[p]; hsv3[2] = col_v[p];
                hsvrgb(hsv3, rgb2);
                col_r[p] = rgb2[0]; col_g[p] = rgb2[1]; col_b[p] = rgb2[2];
                col2Str[0] = '0' + (int)col_s[p]/100;
                col2Str[1] = '0' + ((int)col_s[p]/10)%10;
                col2Str[2] = '0' + (int)col_s[p]%10; 
            } else {
                col_g[p] = mx-25-x;
                col2Str[0] = '0' + (int)col_g[p]/100;
                col2Str[1] = '0' + ((int)col_g[p]/10)%10;
                col2Str[2] = '0' + (int)col_g[p]%10; 
            };
            if (colIn[0]==1) updateColTheme();
            colIn[1] = 2;
        } else if (mx >= 25+x && mx <= 25+255+x && ( (my >= 135+y && my <= 146+y) || colIn[1]==3) ) {
            if (colIn[5]) {
                col_h[p] = (mx-25-x) *1.41137;
                hsv3[0] = col_h[p]; hsv3[1] = col_s[p]; hsv3[2] = col_v[p];
                hsvrgb(hsv3, rgb2);
                col_r[p] = rgb2[0]; col_g[p] = rgb2[1]; col_b[p] = rgb2[2];
                col1Str[0] = '0' + (int)col_h[p]/100;
                col1Str[1] = '0' + ((int)col_h[p]/10)%10;
                col1Str[2] = '0' + (int)col_h[p]%10; 
            } else {
                col_r[p] = mx-25-x;
                col1Str[0] = '0' + (int)col_r[p]/100;
                col1Str[1] = '0' + ((int)col_r[p]/10)%10;
                col1Str[2] = '0' + (int)col_r[p]%10;
            };
            if (colIn[0]==1) updateColTheme();
            colIn[1] = 3;
        };
    };
};

// To draw tangent line if it is chosen as tracer
void drawTangent(double x0, double y0, int j) {
    double y1, y2; int x;
    
    for(x=0; x<iScreenWidth; x++) {
        if (x==0) y1 = y0 + d_dxSinusoid(x0, j) * (x - x0);
        y2 = y0 + d_dxSinusoid(x0, j) * (x+1 - x0);
        iLine (x, y1, x+1, y2);
        iLine (x, y1+1, x+1, y2+1);
        iLine (x, y1-1, x+1, y2-1);
        y1 = y2;
    };
};

// To draw the grid of the coorrdinate system
void drawGrid (int mode) {
    double i; char str[20];
    
    unit = pow(10, round(log10(iScreenWidth/zoomFac[0])-1));
    if (iScreenWidth/zoomFac[0]/unit < 6) unit*=0.5;
    if (iScreenWidth/zoomFac[0]/unit > 21) unit*=2;

    for (i=origin[0]; i<iScreenWidth; i+=unit*zoomFac[0]) {
        if (mode==0) {
            iSetColor (col_r[100], col_g[100], col_b[100]);
            iLine (i, iScreenHeight, i, 0);
        } else if (mode==1) {
            sprintf(str, "%g", (i-origin[0])/zoomFac[0]);
            if (!highlightSum) iSetColor (col_r[102], col_g[102], col_b[102]);
            iText(i+12, origin[1]+12, str);
        };
    };
    for (i=origin[0]-unit*zoomFac[0]; i>0; i-=unit*zoomFac[0]) {
        if (mode==0) {
            iSetColor (col_r[100], col_g[100], col_b[100]);
            iLine (i, iScreenHeight, i, 0);
        } else if (mode==1) {
            sprintf(str, "%g", (i-origin[0])/zoomFac[0]);
            if (!highlightSum) iSetColor (col_r[102], col_g[102], col_b[102]);
            iText(i+12, origin[1]+12, str);
        };
    };
    for (i=origin[1]; i<iScreenHeight; i+=unit*zoomFac[1]) {
        if (mode==0) {
            iSetColor (col_r[100], col_g[100], col_b[100]);
            iLine (0, i, iScreenWidth, i);
        } else if (mode==1) {
            sprintf(str, "%g", (i-origin[1])/zoomFac[1]);
            if (!highlightSum) iSetColor (col_r[102], col_g[102], col_b[102]);
            iText(origin[0]+12, i+12, str);
        };
    };
    for (i=origin[1]-unit*zoomFac[1]; i>0; i-=unit*zoomFac[1]) {
        if (mode==0) {
            iSetColor (col_r[100], col_g[100], col_b[100]);
            iLine (0, i, iScreenWidth, i);
        } else if (mode==1) {
            sprintf(str, "%g", (i-origin[1])/zoomFac[1]);
            if (!highlightSum) iSetColor (col_r[102], col_g[102], col_b[102]);
            iText(origin[0]+12, i+12, str);
        };
    };
};

// To apply properties obtained from Fourier series calculations to the curves, and to execute the calculations
void apply_fourier (int type_G, double n_P, int deg_P, double c_P[], double b_G, double a, double p) {

    int i, Lchange=0; double an, bn, f_amp, f_phase, L1, L2, I, J;
    n = 98; L1 = -10; L2 = 10; G_fr = 1; G_amp = 1;
    sprintf (G_ampStr, "%lg", G_amp);
    sprintf (G_frStr, "%lg", G_fr);

    #define F_      f(
    #define _F      , 0, type_G, n_P, deg_P, c_P, b_G, a, p, 0, 1)
    
    if ( !isnan(F_ 0.01 _F) && F_ 0.01 _F < 80 && F_ 0.01 _F > -80 ) {
        for(I=0.01; I<10; I+=0.01) {
            if ( isnan(F_ I _F) || F_ I _F > 80 || F_ I _F < -80 ) {
                L2 = I-0.01;
                break;
            } else L2 = 10;
        };
        if ( isnan(F_ 0 _F) || F_ 0 _F > 80 || F_ 0 _F < -80 ) L1 = 0.01;
        else if ( !isnan(F_ -0.01 _F) && F_ -0.01 _F < 80 && F_ -0.01 _F > -80 ) {
            for(I=-0.01; I>-10; I-=0.01) {
                if ( isnan(F_ I _F) || F_ I _F > 80 || F_ I _F < -80 ) {
                    L1 = I+0.01;
                    break;
                } else L1 = -10;
            };
        } else L1 = 0;
    } else if ( !isnan(F_ -0.01 _F)  && F_ -0.01 _F < 80 && F_ -0.01 _F > -80  ) {
        L2 = -0.01;
        for(I=-0.01; I>-10; I-=0.01) {
            if ( isnan(F_ I _F) || F_ I _F > 80 || F_ I _F < -80 ) {
                L1 = I+0.01;
                break;
            } else L1 = -10;
        };
    } else {
        for(I=0; I<20000; I+=0.01) {
            if ( !isnan(F_ I _F) && F_ I _F < 80 && F_ I _F > -80 ) {
                L1 = I+0.01;
                for(J=L1; J<L1+20; J+=0.01) {
                    if ( isnan(F_ J _F) || F_ J _F > 80 || F_ J _F < -80 ) {
                        L2 = J - 0.01;
                        Lchange = 1;
                        break;
                    } else {
                        L2 = L1+20;
                        Lchange = 1;
                    };
                };
                break;
            };
        };
        if (!Lchange) {
            for(I=0; I>-20000; I-=0.01) {
                if ( !isnan(F_ I _F) && F_ I _F < 80 && F_ I _F > -80  ) {
                    L2 = I - 0.01;
                    for(J=L2; J>L2-20; J-=0.01) {
                        if ( isnan(F_ J _F) || F_ J _F > 80 || F_ J _F < -80 ) {
                            L1 = J + 0.01;
                            break;
                        } else L1 = L2-20;
                    };
                    break;
                } else L2 = 0;
            };
        };
    };
    an = integral(L1, L2, 0, type_G, n_P, deg_P, c_P, b_G, a, p, 0)/(L2-L1);
    elev[5] = an;

    for (i=1; i<=98; i++) {
        an = 2 * integral(L1, L2, 0, type_G, n_P, deg_P, c_P, b_G, a, p, 1.0*i) / (L2-L1);       
        bn = 2 * integral(L1, L2, 1, type_G, n_P, deg_P, c_P, b_G, a, p, 1.0*i) / (L2-L1);
        f_amp = sqrt(an*an+bn*bn);
        f_phase = atan2(bn,an);
        amp[i] = f_amp;
        fr[i] = (2*PI/(L2-L1)) * i*180/PI;
        phase[i] = -f_phase*180/PI;
        sincos[i] = 90;
        tracer_hid[i] = 1;
        tracer_x[i] = tracer_phase[i] = 27;
        path_hid[i] = 0;
        if (i!=5) elev[i] = 0;
        col_r[i] = col_r[100];  col_g[i] = col_g[100];  col_b[i] = col_b[100];
        tweakIn[4] = i;
        tweakDispInit();
    };
    highlightSum = 1;
    updateColTheme();
    zoomFac[0] = zoomFac[1] = 63.80912603;
    origin[0]=iScreenWidth/2; origin[1]=190;
};

// ____________________THE MAIN iDraw() FUNCTION ________________
void iDraw() {

    int i, j;
    iClear();
    iSetColor (col_r[101], col_g[101], col_b[101]);
    iFilledRectangle(0, 0, iScreenWidth, iScreenHeight);
    
    if (!path_hid[0]) {
        drawGrid(0);        // Draw the grid of coordinate axes
        // Draw the arrrows at end of coordinate axxes, if the other axis isn't to close to screen edge
        iSetColor(col_r[102], col_g[102], col_b[102]);
        double x_x[3]={0,30,30}, y_x[3]={origin[1], origin[1]+10, origin[1]-10}, yy[3]={iScreenHeight-30, iScreenHeight, iScreenHeight-30},
        xx[3]={iScreenWidth,iScreenWidth-30,iScreenWidth-30}, x_y[3]={origin[0]-10,origin[0],origin[0]+10}, y_y[3]={30,0,30};
        if (origin[0]>100) iFilledPolygon (x_x, y_x, 3); 
        if (origin[0]<iScreenWidth-100) iFilledPolygon (xx, y_x, 3); 
        if (origin[1]<iScreenHeight-100) iFilledPolygon (x_y, yy, 3); 
        if (origin[1]>100) iFilledPolygon (x_y, y_y, 3);
        // Draw bold coordinate axes
        iLine (0, origin[1]-1, iScreenWidth, origin[1]-1);
        iLine (0, origin[1], iScreenWidth, origin[1]);
        iLine (0, origin[1]+1, iScreenWidth, origin[1]+1);
        iLine (origin[0]-1, iScreenHeight, origin[0]-1, 0);
        iLine (origin[0], iScreenHeight, origin[0], 0);
        iLine (origin[0]+1, iScreenHeight, origin[0]+1, 0);

        #define AMP_FR_RATIO  !((i+pathAnimState)%5) && fabs(amp[j]*fr[j]*G_amp*G_fr)<515.05  ||  blinkAnimState && fabs(amp[j]*fr[j]*G_amp*G_fr)>=515.05

        for(i=0; i<iScreenWidth; i++) {
            
            if (dot) {                          // draw curves and some in dot plotting mode
                for (j=1; j<=n; j++) {
                    iSetColor (col_r[j], col_g[j], col_b[j]);
                    y_i[j] = sinusoid(i, j);
                    if (!path_hid[j]) {
                        if ( tweakIn[4]==j && tweakIn[0]!=0 ) {
                            //Make curve bold and display chain animation when open in edit menu
                            if ( AMP_FR_RATIO ) {
                                iPoint (i-2, y_i[j]);
                                iPoint (i-1, y_i[j]);
                                iPoint (i+1, y_i[j]);
                                iPoint (i+2, y_i[j]);
                                iPoint (i, y_i[j]);
                                iPoint (i, y_i[j]-1);
                                iPoint (i, y_i[j]+1);
                                iPoint (i, y_i[j]-2);
                                iPoint (i, y_i[j]+2);
                            };
                        } else {
                            iPoint (i, y_i[j]);
                        };
                        if ( drag==j ) {
                            // Make curves bold when selected/dragging
                            iPoint (i-2, y_i[j]);
                            iPoint (i-1, y_i[j]);
                            iPoint (i+1, y_i[j]);
                            iPoint (i+2, y_i[j]);
                            iPoint (i, y_i[j]-1);
                            iPoint (i, y_i[j]+1);
                            iPoint (i, y_i[j]-2);
                            iPoint (i, y_i[j]+2);
                        };
                    };
                };
            } else {                                      // draw curves and some in continuous line mode
                for (j=1; j<=n; j++) {
                    iSetColor (col_r[j], col_g[j], col_b[j]);
                    if(i==0) y_i[j] = sinusoid(i, j);
                    y_ip1[j] = sinusoid(i,j);
                    if (!path_hid[j]) {
                        if ( tweakIn[4]==j && tweakIn[0]!=0 ) {
                            //Make curve bold and display chain animation when open in edit menu
                            if ( AMP_FR_RATIO ) {
                                iLine (i-2, y_i[j], i-1, y_ip1[j]);
                                iLine (i-1, y_i[j], i, y_ip1[j]);
                                iLine (i+1, y_i[j], i+2, y_ip1[j]);
                                iLine (i+2, y_i[j], i+3, y_ip1[j]);
                                iLine (i, y_i[j], i+1, y_ip1[j]);           //----------
                                iLine (i, y_i[j]-1, i+1, y_ip1[j]-1);
                                iLine (i, y_i[j]+1, i+1, y_ip1[j]+1);
                                iLine (i, y_i[j]-2, i+1, y_ip1[j]-2);
                                iLine (i, y_i[j]+2, i+1, y_ip1[j]+2);
                            };
                        } else {
                            iLine (i, y_i[j], i+1, y_ip1[j]);
                        };
                        if ( (drag==j) ) {
                            // Make curves bold when selected/dragging
                            iLine (i-2, y_i[j], i-1, y_ip1[j]);
                            iLine (i-1, y_i[j], i, y_ip1[j]);
                            iLine (i+1, y_i[j], i+2, y_ip1[j]);
                            iLine (i+2, y_i[j], i+3, y_ip1[j]);           //....
                            iLine (i, y_i[j]-1, i+1, y_ip1[j]-1);
                            iLine (i, y_i[j]+1, i+1, y_ip1[j]+1);
                            iLine (i, y_i[j]-2, i+1, y_ip1[j]-2);
                            iLine (i, y_i[j]+2, i+1, y_ip1[j]+2);
                        };
                    };
                };
            };

            iSetColor (col_r[0], col_g[0], col_b[0]);
            if (!dot) {                 // Draw summaiton curves in dot mode
                iLine ( i, arrsum(y_i)-(origin[1])*(n-1), i+1, arrsum(y_ip1)-(origin[1])*(n-1) );
                iLine ( i+1, arrsum(y_i)-(origin[1])*(n-1), i+2, arrsum(y_ip1)-(origin[1])*(n-1) );
                iLine ( i, arrsum(y_i)-(origin[1])*(n-1)-1, i+1, arrsum(y_ip1)-(origin[1])*(n-1)-1 );
                if (highlightSum) {         //Make highlighted if highlight mode slected
                    iLine ( i-2, arrsum(y_i)-(origin[1])*(n-1), i-1, arrsum(y_ip1)-(origin[1])*(n-1) );
                    iLine ( i-1, arrsum(y_i)-(origin[1])*(n-1), i, arrsum(y_ip1)-(origin[1])*(n-1) );
                    iLine ( i+1, arrsum(y_i)-(origin[1])*(n-1), i+2, arrsum(y_ip1)-(origin[1])*(n-1) );
                    iLine ( i+2, arrsum(y_i)-(origin[1])*(n-1), i+3, arrsum(y_ip1)-(origin[1])*(n-1) );
                    iLine ( i, arrsum(y_i)-(origin[1])*(n-1)-2, i+1, arrsum(y_ip1)-(origin[1])*(n-1)-2 );
                    iLine ( i, arrsum(y_i)-(origin[1])*(n-1)-1, i+1, arrsum(y_ip1)-(origin[1])*(n-1)-1 );
                    iLine ( i, arrsum(y_i)-(origin[1])*(n-1)+1, i+1, arrsum(y_ip1)-(origin[1])*(n-1)+1 );
                    iLine ( i, arrsum(y_i)-(origin[1])*(n-1)+2, i+1, arrsum(y_ip1)-(origin[1])*(n-1)+2 );
                };
                for(j=1; j<=n; j++) y_i[j]=y_ip1[j];
            } else {            //Draw summation curve in line mode
                iPoint ( i, arrsum(y_i)-(origin[1])*(n-1) );
                if (highlightSum) {         //Make highlighted if highlight mode slected
                    iPoint ( i-2, arrsum(y_i)-(origin[1])*(n-1) );
                    iPoint ( i-1, arrsum(y_i)-(origin[1])*(n-1) );
                    iPoint ( i+1, arrsum(y_i)-(origin[1])*(n-1) );
                    iPoint ( i+2, arrsum(y_i)-(origin[1])*(n-1) );
                    iPoint ( i, arrsum(y_i)-(origin[1])*(n-1)-2 );
                    iPoint ( i, arrsum(y_i)-(origin[1])*(n-1)-1 );
                    iPoint ( i, arrsum(y_i)-(origin[1])*(n-1)+1 );
                    iPoint ( i, arrsum(y_i)-(origin[1])*(n-1)+2 );
                };
            };

            // To show the actual graph of the function along with the fourier series plot, when in fourier plot mode
            if (f_mode==2) {
                iSetColor (249, 12, 14);
                if(i==0) y_i[0] = ( f( (i-origin[0])/zoomFac[0], 0, type_G, n_P, deg_P, c_arr, b_G, a, p, 0, 1) )*zoomFac[1] + origin[1];
                y_ip1[0] = ( f( (i+1-origin[0])/zoomFac[0], 0, type_G, n_P, deg_P, c_arr, b_G, a, p, 0, 1) )*zoomFac[1] + origin[1];
                if (!isnan(y_i[0]) && !isnan(y_ip1[0]) ) {
                    iLine (i, y_i[0], i+1, y_ip1[0] );
                };
                y_i[0] = y_ip1[0];
            };

        };        

    };

    // Draw tracers of curves 
    if (!tracer_hid[0]) {
        for (j=1;j<=n;j++) {
            if (!tracer_hid[j]) {
                iSetColor (col_r[j], col_g[j], col_b[j]);
                if (tracer_img[j]==0) iFilledCircle ( tracer_x[j], sinusoid(tracer_x[j],j), 6 );
                else if (tracer_img[j]==1) {
                    iRotate(tracer_x[j], sinusoid(tracer_x[j],j), 180/PI*atan(d_dxSinusoid(tracer_x[j],j)));
                    iFilledRectangle ( tracer_x[j]-8, sinusoid(tracer_x[j],j)-8, 16, 16);
                    iUnRotate();
                } else if (tracer_img[j]==9) {
                    iFilledCircle ( tracer_x[j], sinusoid(tracer_x[j],j), 5 );
                    drawTangent(tracer_x[j], sinusoid(tracer_x[j],j), j);
                } else if (tracer_img[j]==10) {
                    iRotate(tracer_x[j], sinusoid(tracer_x[j],j), 180/PI*atan(d_dxSinusoid(tracer_x[j],j)));
                    drawPacman (blinkAnimState, tracer_x[j], sinusoid(tracer_x[j], j), j );
                    iUnRotate();
                } else iShowBMP2 (tracer_x[j]-25, sinusoid(tracer_x[j],j)-25, imgNames[tracer_img[j]], MATTE);
            };
        };

        if (n>1) {
            iSetColor (col_r[0], col_g[0], col_b[0]);
            sinosum = origin[1];
            for (i=1; i<=n; i++) sinosum += (sinusoid(tracer_x[0], i) - origin[1]);
            iFilledCircle( tracer_x[0], sinosum , (9+3*highlightSum) );
            sinosum = origin[1];
        }; 
    };

    // Draw graduations on coordinate axes, over the curves
    drawGrid(1);

    drawButtons();
    drawGFacs();

    if (colIn[0]) drawColMenu (colIn[2], colIn[3]); 
    if (tweakIn[0]) drawTweakMenu (tweakIn[2], tweakIn[3]);
    if (delIn[0]) drawDelMenu (delIn[1], delIn[2]);
    if (f_mode==1) drawFourierHelp();
    if (help==1) iShowBMP (0,0,"help.bmp");

};


void iMouseMove(int mx, int my) {
    if (MB == GLUT_LEFT_BUTTON) {
        // Respond to dragging curves, when the correct value of MB is set by iMouse()
        if (h1>origin[0]) fr[drag] = 360/( 360/fr2[drag] + sgn(fr2[drag])*(1.0*mx-h1)/zoomFac[0] );
        else fr[drag] = 360/( 360/fr2[drag] - sgn(fr2[drag])*(1.0*mx-h1)/zoomFac[0] );
        if (v1>elev[drag]+origin[1]) amp[drag] = amp2[drag] + sgn(amp2[drag])*(1.0*my-v1)/zoomFac[1];
        else amp[drag] = amp2[drag] - sgn(amp2[drag])*(1.0*my-v1)/zoomFac[1];
    } else if (MB == GLUT_RIGHT_BUTTON) {
        // Respond to dragging curves, when the correct value of MB is set by iMouse()
        phase[drag] = phase2[drag] - sgn(fr[drag])*(1.0*mx-h1)/zoomFac[0];
        elev[drag] = elev2[drag] + (1.0*my-v1)/zoomFac[1];
    } else if (MB == GLUT_MIDDLE_BUTTON) {
        // Respond to panning with MMB, when the correct value of MB is set by iMouse()
        if (my-v1 > 0) if (origin[1]<1400000*zoomFac[1]) origin[1] = origin2[1] + (my-v1);
        if (my-v1 < 0) if (origin[1]>-175000*zoomFac[1]) origin[1] = origin2[1] + (my-v1);
        if (mx-h1 > 0) if (origin[0]<2100000*zoomFac[0]) origin[0] = origin2[0] + (mx-h1);
        if (mx-h1 < 0) if (origin[0]>-175000*zoomFac[0]) origin[0] = origin2[0] + (mx-h1);
    };
    // Update values of colour or curve properties if sliders are moved in add/edit or colour input menus
    if (colIn[0]) colValue (mx, my, colIn[2], colIn[3], colIn[4], GLUT_DOWN);
    if (tweakIn[0]) tweakValue (mx, my, tweakIn[2], tweakIn[3], GLUT_DOWN);
};


void iMouse(int button, int state, int mx, int my) {

    if (f_mode!=1 && help==0) {                 // If the Fourier series menu of help screen aren't covering the screen
        //Defining the areas for UI elements on screen
        #define T1 (mx>=138 && mx<=478 && my>=iScreenHeight-296 && my<=iScreenHeight-24)
        #define T2 (mx>=78 && mx<=418 && my>=iScreenHeight-296 && my<=iScreenHeight-24)
        #define C1 (mx>=iScreenWidth-536 && mx<=iScreenWidth-196 && my>=iScreenHeight-184 && my<=iScreenHeight-24)
        #define C2 (mx>=424 && mx<=764 && my>=iScreenHeight-248 && my<=iScreenHeight-88)
        #define C2P (mx>=484 && mx<=824 && my>=iScreenHeight-248 && my<=iScreenHeight-88)
        #define VT1 (tweakIn[0]==1)
        #define VT2 (tweakIn[0]==2)
        #define VC1 (colIn[0]==1)
        #define VC2 (colIn[0]==2)
        #define D1 (delIn[0] && mx>=324 && mx<=524 && my>=iScreenHeight-160 && my<=iScreenHeight-65)
        #define B1 (mx>=24 && mx<=372 && my>=iScreenHeight-60 && my<=iScreenHeight-24)
        #define B2 (mx>=iScreenWidth-192 && mx<=iScreenWidth-24 && my>=iScreenHeight-60 && my<=iScreenHeight-24)
        #define B3 (mx>=24 && mx<=132 && my>=24 && my<=60)
        #define MOUSE_SELECT (fabs(my-sinusoid(mx, i)) <= 10  ||  fabs(fmod(mx,90)-fmod(inverseSinusoid(my, i),90)) <= 10)
        #define NO_DRAG (drag2 == 0)

        int i; MB = button;
        if (button == GLUT_LEFT_BUTTON){
            if(state == GLUT_DOWN) {
                // Change values of properties when curves are dragged, but the mouse pointer is NOT in the area of any UI element
                if (!drag && !(  D1||B1||B2||B3||(VT1&&T1)||(VT2&&T2)||(VC1&&C1)||(VC2&&((VT1&&C2P)||(VT2&&C2)))  ) ) {
                    for (i=1; i<=n; i++) {
                        if( MOUSE_SELECT && path_hid[i]==0) {
                            // if the point almost matches a point on the set of solutions of any curve, fix the drag index to that curves index
                            drag = i;
                            break;
                        } else drag =0;
                    };
                    fr2[i] = fr[i]; amp2[i] = amp[i];
                    h1 = mx; v1 = my;
                    drag2 = drag;
                };

                if (colIn[0]) colValue(mx,my,colIn[2], colIn[3], colIn[4], state);
                if (tweakIn[0]) tweakValue (mx, my, tweakIn[2], tweakIn[3], state);
                
            } else if(state == GLUT_UP) { 
                // Only if no curves are being dragged, set drag2 to 0
                // because if a curve starts reponding to dragging, those are set to the coordinates of the mouse pointer
                if (!drag) drag2 = 0;
                // set the drag index to zero as no curve is being drargged
                drag = 0;
                // If colour or edit/add menus are open, then update values of quantities that were being changed by dagging sliders
                if (colIn[0]) {
                    colIn[1] = 0;
                    if (colIn[0]) colValue(mx,my,colIn[2], colIn[3], colIn[4], state);
                };
                if (tweakIn[0]) {
                    tweakIn[1] = 0;
                    tweakValue (mx, my, tweakIn[2], tweakIn[3], state);
                    if ( tweakIn[0]==1 && !(  D1||B1||B2||B3||(VT1&&T1)||(VT2&&T2)||(VC1&&C1)||(VC2&&((VT1&&C2P)||(VT2&&C2)))  ) ) {
                        for (i=1; i<=n; i++) {
                            if( MOUSE_SELECT && path_hid[i]==0 && drag2==i) {
                                // if the point almost matches a point on the set of solutions of any curve, fix the drag index to that curves index
                                tweakIn[4] = i;
                                break;
                            };
                        };
                        tweakDispInit();
                        if (colIn[0] == 2) {
                            colIn[4] = tweakIn[4];
                            colDispInit(colIn[4]);
                        };
                    };
                };        
                //__________Responses trigerred by releasing LMB in areas of UI buttons, but not if a curve was being dragged        
                if (NO_DRAG && mx >= iScreenWidth-72 && mx<=iScreenWidth-24 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    if (dot) dot = 0;
                    else dot = 1;
                } else if (NO_DRAG && mx >= iScreenWidth-132 && mx<=iScreenWidth-84 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    if (cont_mode) cont_mode = 0;
                    else cont_mode = 1;
                } else if (NO_DRAG && mx >= iScreenWidth-192 && mx<=iScreenWidth-144 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    colIn[2] = iScreenWidth-536; colIn[3] = iScreenHeight-184; 
                    if (colIn[0]) {     
                        if (colIn[0]==2) colIn[0]=0;
                        else {               
                            colIn[0] = 0;
                        };
                    } else {          
                        colIn[4] = 101;
                        colDispInit(101);
                        colIn[0] = 1; 
                    };
                } else if (NO_DRAG && mx >= 24 && mx<=72 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    delIn[0] = 0;
                    if (tweakIn[0]) { 
                        if (colIn[0]) colIn[0]=0;                   
                        if (tweakIn[0]==2) {
                            tweakIn[0] = 0;
                            deleteCurve(n);
                            tweakDispInit();
                        } else tweakIn[0] = 0;
                    } else {          
                        if (n<98) {
                            tweakIn[2] = 78; tweakIn[3] = iScreenHeight-224;
                            tweakIn[0] = 2;
                            addCurve();
                            tweakDispInit();
                        };      
                    };
                } else if (NO_DRAG && tweakIn[0]!=2 && mx >= 84 && mx<=132 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    delIn[0] = 0;
                    if (tweakIn[0]) {                    
                        tweakIn[0] = 0;
                        if (colIn[0]) colIn[0]=0;
                    } else if (n>0) {
                        tweakIn[2] = 138; tweakIn[3] = iScreenHeight-224;
                        tweakIn[0] = 1;
                        tweakDispInit();
                    };
                } else if (NO_DRAG && !tweakIn[0] && mx >= 144 && mx<=192 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    int k;
                    for (k=0; k<=n; k++) {
                        tracer_x[k] = tracer_phase[k];
                        dir[k] = 1;
                    };
                } else if (NO_DRAG && !tweakIn[0] && mx >= 204 && mx<=252 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    if (highlightSum) highlightSum = 0;
                    else highlightSum = 1;
                    updateColTheme();
                } else if (NO_DRAG && mx >= 24 && mx<=72&& my>= 24 && my<= 60) {
                    origin[0] = iScreenWidth/2;
                    origin[1] = iScreenHeight/2;
                } else if (NO_DRAG && mx >= 84 && mx<=132 && my>= 24 && my<= 60) {
                    textIn = 0;
                    help =1;
                } else if (NO_DRAG && !tweakIn[0] && mx >= 324 && mx<=372 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    if (!delIn[0]) {
                        delIn[1] = 324; delIn[2]=iScreenHeight-60;
                        delIn[0] = 1; 
                    } else delIn[0] = 0;
                } else if (NO_DRAG && delIn[0]==1 && mx >= 334 && mx<=419 && my>= iScreenHeight-150 && my<= iScreenHeight-114) {
                    delIn[0] = 0; n = 0;
                    if (f_mode==2) f_mode = 0;
                    tweakDispInit();
                } else if (NO_DRAG && delIn[0]==1 && mx >= 429 && mx<=514 && my>= iScreenHeight-150 && my<= iScreenHeight-114) {
                    delIn[0] = 0; 
                } else if (NO_DRAG && !tweakIn[0] && mx >= 264 && mx<=312 && my>= iScreenHeight-60 && my<= iScreenHeight-24) {
                    textIn = 0;
                    fourierDispInit();
                    if (f_mode!=1) f_mode = 1;
                } else if (NO_DRAG && mx >= iScreenWidth-92 && mx<=iScreenWidth-78 && my>= 10 && my<= 20) {
                    if (abs(dx)>1.5) (dx>0)? dx-- : dx++;
                    else if (abs(dx)>0.2) (dx>0)? dx-=0.1 : dx+=0.1;
                    snprintf(dxStr, 18, "%.2lf", dx);
                } else if (NO_DRAG && mx >= iScreenWidth-46 && mx<=iScreenWidth-32 && my>= 10 && my<= 20) {
                    if(abs(dx)<50) (dx>0)? dx++ : dx--;
                    snprintf(dxStr, 18, "%.2lf", dx);
                };                
            };
        } else if (button == GLUT_RIGHT_BUTTON){
            if (state == GLUT_DOWN) {
                if (!drag && !(  D1||B1||B2||B3||(VT1&&T1)||(VT2&&T2)||(VC1&&C1)||(VC2&&((VT1&&C2P)||(VT2&&C2)))  ) ) {
                    for (i=1; i<=n; i++) {
                        if( MOUSE_SELECT && path_hid[i]==0) {
                            drag = i;
                            break;
                        } else drag =0;
                    };
                    phase2[i] = phase[i]; elev2[i] = elev[i];
                    h1 = mx; v1 = my;
                    drag2 = drag;
                };
            } else if(state == GLUT_UP) { 
                if (!drag) drag2 = 0;
                drag = 0;
                if (NO_DRAG && mx >= iScreenWidth-192 && mx<=iScreenWidth-144 && my>= iScreenHeight-60 && my<= iScreenHeight-24 && colIn[0]!=2) {
                    double rgb1[6], hsv1[6];
                    col_r[101] = 28;
                    col_g[101] = 33;
                    col_b[101] = 47;
                    rgb1[0] = col_r[101]; rgb1[1] = col_g[101]; rgb1[2] = col_b[101];
                    rgbhsv(rgb1, hsv1);
                    col_h[101]=hsv1[0]; col_s[101]=hsv1[1]; col_v[101]=hsv1[2];
                    colDispInit(101);
                    updateColTheme();
                };
                if (tweakIn[0]) {
                    tweakIn[1] = 0;
                    if (tweakIn[0]) tweakValue (mx, my, tweakIn[2], tweakIn[3], state);
                    tweakDispInit();
                };                    
            };
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (state == GLUT_DOWN) {
                origin2[0] = origin[0];
                origin2[1] = origin[1];
                h1 = mx; v1 = my;
            };
        };
    } else if (f_mode==1 && state==GLUT_UP && button==GLUT_LEFT_BUTTON) {
        // Button and textbox responses for clicks in the Fourier series unction input menu
        if (mx >= 280 && mx<=294 && my>= 220 && my<= 244) {
            if (type_G != 20) type_G++;
            else type_G=0;
            fourierDispInit();
        } else if (mx >= 220 && mx<=234 && my>= 220 && my<= 244) {
            if (type_G != 0) type_G--;
            else type_G=20;
            fourierDispInit();
        } else if (mx >= 340 && mx<=394 && my>= 220 && my<= 244) {
            textIn = 6;
            fourierDispInit();
        } else if (mx >= 448 && mx<=502 && my>= 220 && my<= 244) {
            textIn = 7;
            fourierDispInit();
        } else if (mx >= 588 && mx<=642 && my>= 220 && my<= 244) {
            textIn = 8;
            fourierDispInit();
        } else if (mx >= 850 && mx<=904 && my>= 220 && my<= 244) {
            textIn = 9;
            fourierDispInit();
        } else if (mx >= 788 && mx<=1168 && my>= 178 && my<= 202) {
            textIn = 10;
            fourierDispInit();
        } else if (mx >= 1018 && mx<=1072 && my>= 220 && my<= 244) {
            textIn = 11;
            fourierDispInit();
        } else textIn = 0;
        if (mx >= iScreenWidth/2-32 && mx<=iScreenWidth/2+44 && my>= 24 && my<= 60) {
            textIn = 0;            
            apply_fourier(type_G, n_P, deg_P, c_arr, b_G, a, p);
            f_mode = 2;
            textIn = 0;
            fourierDispInit();
        } else if (mx >= iScreenWidth-56 && mx<=iScreenWidth-24 && my>= iScreenHeight-56 && my<= iScreenHeight-24) {
            f_mode = 0;
            textIn = 0;
            fourierDispInit();     
        };
    } else if (help==1 && state==GLUT_UP && button==GLUT_LEFT_BUTTON) {
        // Button and textbox responses for clicks in the help screen
        if (mx>=16 && mx<=45 && my>=iScreenHeight-43 && my<=iScreenHeight-15) {
            help=0;
        };
    };
};


void iKeyboard(unsigned char key) {
    if (f_mode!=1 && help==0) {
        if (key == 'p' || key == 'P') {
            if (!paused[0]) {
                iPauseTimer(0);
                paused[0] = 1;
            } else {
                iResumeTimer(0);
                paused[0] = 0;
            };
        };
        if (key == 'r' || key == 'R') iResumeTimer(0);
        if (key == 's') {
            if (path_hid[0]) path_hid[0] = 0;
            else path_hid[0] = 1;
        };
        if (key == 'S') if (path_hid[0]) path_hid[0] = 0;
        if (key == '-' && !textIn) {
            if (abs(dx)>1.5) (dx>0)? dx-- : dx++;
            else if (abs(dx)>0.2) (dx>0)? dx-=0.1 : dx+=0.1;
            snprintf(dxStr, 18, "%.2lf", dx);
        };
        if ((key == '=' || key =='+') && !textIn) {
            if(abs(dx)<50) (dx>0)? dx++ : dx--;
            snprintf(dxStr, 18, "%.2lf", dx);
        };
        if (key == 'a' && G_amp<6.33) {
            G_amp += 0.02;
            sprintf (G_ampStr, "%lg", G_amp);
        };
        if (key == 'A' && G_amp>0) {
            G_amp -= 0.02;
            sprintf (G_ampStr, "%lg", G_amp);
        };
        if (key == 'f' && G_fr<4) {
            G_fr += 0.01;
            sprintf (G_frStr, "%lg", G_fr);
        };
        if (key == 'F' && G_fr>0) {
            G_fr -= 0.01;
            sprintf (G_frStr, "%lg", G_fr);
        };
        if (key>='0' && key<='9' && !textIn) {
            if (tracer_hid[key-'0']) tracer_hid[key-'0'] = 0;
            else tracer_hid[key-'0'] = 1;
        };
        if (key == 'c') {
            if (cont_mode) cont_mode = 0;
            else cont_mode = 1;
        };
        if (key == 'C') if (cont_mode) cont_mode = 0;
        if (key == 'd') {
            if (dot) dot = 0;
            else dot = 1;
        };
        if (key == 'D') if (dot) dot = 0;

    };

    if (key == 6) {
        textIn = 0; help=0;
        fourierDispInit();
        if (f_mode==1) f_mode = 0;
        else f_mode = 1;
    };
    if (key == 'q') exit(0);

    // ____________Text box keyboard inputs in the edit/add menu_________
    if (textIn==1) {          
        if (key == '\r') {
            textIn = 0;
            sscanf (tweakNoStr, "%d", &tweakIn[4]);
            if (tweakIn[4] < 1) tweakIn[4] = 1;
            if (tweakIn[4] > n) tweakIn[4] = n;
            tweakDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(tweakNoStr) > 0) tweakNoStr[strlen(tweakNoStr)-1] = '\0';
        } else if (key>='0' && key<='9') {
            if (strlen(tweakNoStr) < 3) {                
                tweakNoStr[strlen(tweakNoStr)+1] = '\0';
                tweakNoStr[strlen(tweakNoStr)] =  key;
            };
        };
    } else if (textIn==2) {          
        if (key == '\r') {
            sscanf (ampStr, "%lf", &amp[tweakIn[4]]);
            if (amp[tweakIn[4]]>500000) amp[tweakIn[4]] = 500000;
            else if (amp[tweakIn[4]]<-500000) amp[tweakIn[4]] = -500000;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(ampStr) > 0) ampStr[strlen(ampStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(ampStr) < 8) {
                ampStr[strlen(ampStr)+1] = '\0';                
                ampStr[strlen(ampStr)] =  key;
            };
        };
    } else if (textIn==3) {          
        if (key == '\r') {
            sscanf (frStr, "%lf", &fr[tweakIn[4]]);
            if (fr[tweakIn[4]]>11300) fr[tweakIn[4]] = 11300;
            else if (fr[tweakIn[4]]<-11300) fr[tweakIn[4]] = -11300;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(frStr) > 0) frStr[strlen(frStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(frStr) < 8) {
                frStr[strlen(frStr)+1] = '\0';
                frStr[strlen(frStr)] =  key;                
            };
        };
    } else if (textIn==4) {          
        if (key == '\r') {
            sscanf (elevStr, "%lf", &elev[tweakIn[4]]);
            if (elev[tweakIn[4]]>200000) elev[tweakIn[4]] = 200000;
            else if (elev[tweakIn[4]]<-200000) elev[tweakIn[4]] = -200000;
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(elevStr) > 0) elevStr[strlen(elevStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(elevStr) < 8) {
                elevStr[strlen(elevStr)+1] = '\0';
                elevStr[strlen(elevStr)] =  key;                
            };
        };
    } else if (textIn==5) {          
        if (key == '\r') {
            sscanf (phaseStr, "%lf", &phase[tweakIn[4]]);
            phase[tweakIn[4]] = fmod(phase[tweakIn[4]], 360);
            if (f_mode==2) f_mode = 0;
            tweakDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(phaseStr) > 0) phaseStr[strlen(phaseStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(phaseStr) < 8) {
                phaseStr[strlen(phaseStr)+1] = '\0';
                phaseStr[strlen(phaseStr)] =  key;
            };
        };

    // ________________Text box keyboard inputs in the Fourier series menu_______________
    } else if (textIn==6) {          
        if (key == '\r') {
            sscanf (aStr, "%lf", &a);
            if (a>10e100) a=10e100;
            if (a<-10e100) a=-10e100;
            fourierDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(aStr) > 0) aStr[strlen(aStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(aStr) < 10) {
                aStr[strlen(aStr)+1] = '\0';
                aStr[strlen(aStr)] =  key;
            };
        };
    } else if (textIn==7) {          
        if (key == '\r') {
            sscanf (pStr, "%lf", &p);
            if (p>20) p=20;
            if (p<-20) p=-20;
            fourierDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(pStr) > 0) pStr[strlen(pStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(pStr) < 10) {
                pStr[strlen(pStr)+1] = '\0';
                pStr[strlen(pStr)] =  key;
            };
        };
    } else if (textIn==8) {          
        if (key == '\r') {
            sscanf (n_PStr, "%lf", &n_P);
            if (n_P>20) n_P=20;
            if (n_P<-20) n_P=-20;
            fourierDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(n_PStr) > 0) n_PStr[strlen(n_PStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(n_PStr) < 8) {
                n_PStr[strlen(n_PStr)+1] = '\0';
                n_PStr[strlen(n_PStr)] =  key;
            };
        };
    } else if (textIn==9) {          
        if (key == '\r') {
            sscanf (deg_PStr, "%d", &deg_P);
            if (deg_P>12) deg_P=12;
            if (deg_P<0) deg_P=0;
            fourierDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(deg_PStr) > 0) deg_PStr[strlen(deg_PStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(deg_PStr) < 8) {
                deg_PStr[strlen(deg_PStr)+1] = '\0';
                deg_PStr[strlen(deg_PStr)] =  key;
            };
        };
    } else if (textIn==10) {  
        char *temp; int count=0, i=0, j, k;    
        if (key == '\r') {
            sscanf (c_arrStr, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &c_arr[0],&c_arr[1],&c_arr[2],&c_arr[3],&c_arr[4],&c_arr[5],&c_arr[6],&c_arr[7],&c_arr[8],&c_arr[9],&c_arr[10],&c_arr[11],&c_arr[12]);
            fourierDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(c_arrStr) > 0) c_arrStr[strlen(c_arrStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.' || key==' ') {
            if (strlen(c_arrStr) < 128) {
                c_arrStr[strlen(c_arrStr)+1] = '\0';
                c_arrStr[strlen(c_arrStr)] =  key;
            };
        };
    } else if (textIn==11) {          
        if (key == '\r') {
            sscanf (b_GStr, "%lf", &b_G);
            if (b_G>10) b_G=10;
            if (b_G<-10) b_G=-10;
            if (type_G==19 && b_G<=0) b_G=1;
            fourierDispInit();
            textIn = 0;
        } else if (key == '\b') {
            if (strlen(b_GStr) > 0) b_GStr[strlen(b_GStr)-1] = '\0';
        } else if (key=='-' || key=='+' || key=='e' || key=='E' || key>='0' && key<='9' || key=='.') {
            if (strlen(b_GStr) < 8) {
                b_GStr[strlen(b_GStr)+1] = '\0';
                b_GStr[strlen(b_GStr)] =  key;
            };
        };
    };
};


void iSpecialKeyboard(unsigned char key) {
    
    if (key == GLUT_KEY_END) exit(0);

    //__________For pan and zoom with keyboard__________

    if (f_mode!=1 && help==0) {
        if (key == GLUT_KEY_PAGE_UP) {
            if (zoomFac[0]<2000 && zoomFac[1]<2000) {
                if (zoomFac[0]>=1 || zoomFac[0]>=1) {
                    zoomFac[0] *=1.08;
                    zoomFac[1] *=1.08;
                    origin[0] -= (0.08*(iScreenWidth/2-origin[0]));
                    origin[1] -= (0.08*(iScreenHeight/2-origin[1]));
                } else {
                    zoomFac[0] *=1.05;
                    zoomFac[1] *=1.05;
                    origin[0] -= (0.05*(iScreenWidth/2-origin[0]));
                    origin[1] -= (0.05*(iScreenHeight/2-origin[1]));
                };
            };
        };
        if (key == GLUT_KEY_PAGE_DOWN) {
            if (zoomFac[0]>0.0005 && zoomFac[1]>0.0005) {
                if (zoomFac[0]>1 || zoomFac[0]>1) {
                    zoomFac[0] *= 0.95;
                    zoomFac[1] *= 0.95;
                    origin[0] += (0.05*(iScreenWidth/2-origin[0]));
                    origin[1] += (0.05*(iScreenHeight/2-origin[1]));
                } else {
                    zoomFac[0] *= 0.98;
                    zoomFac[1] *= 0.98;
                    origin[0] += (0.02*(iScreenWidth/2-origin[0]));
                    origin[1] += (0.02*(iScreenHeight/2-origin[1]));
                };
            };
        };
        if (key == GLUT_KEY_UP) if (origin[1]>-175000*zoomFac[1]) origin[1] -= 5;
        if (key == GLUT_KEY_DOWN) if (origin[1]<1400000*zoomFac[1]) origin[1] += 5;
        if (key == GLUT_KEY_RIGHT) if (origin[0]>-175000*zoomFac[0]) origin[0] -= 5;
        if (key == GLUT_KEY_LEFT) if (origin[0]<2100000*zoomFac[0]) origin[0] += 5;
    };

    if (key == GLUT_KEY_F1) {
        textIn=0;
        f_mode=0; tweakIn[0]=0; colIn[0]=0;
        if (help==0) help = 1;
        else help = 0;
    };
};

// To render and animate the curve tracers
void trace() {
    int i;
    for (i=0; i<=n; i++) tracer_x[i] += (dx*dir[i]);
    if (cont_mode) {
        for (i=0; i<=n; i++) {            
            if (tracer_x[i]>iScreenWidth) tracer_x[i]=0;
            else if (tracer_x[i]<0) tracer_x[i]=iScreenWidth;
        };
    } else {
        for (i=0; i<=n; i++) {
            if (tracer_img[i]>=2 && tracer_img[i]<=8) {
                if(tracer_x[i]>=iScreenWidth-26 || tracer_x[i]<=26) dir[i] = -dir[i];
            } else {
                if(tracer_x[i]>iScreenWidth || tracer_x[i]<0) dir[i] = -dir[i];
            };
        };
    };
};

// For the chain animation of the curve that is selected in add/edit menu
void pathAnim() {
    if (pathAnimState==5) pathAnimState=0;
    pathAnimState++;
};

// For any blinking animation conditioning, and for the opening/closing of Pacman's mouth 
void blinkAnim() {
    if (blinkAnimState) blinkAnimState = 0;
    else blinkAnimState = 1;
};


int main() {

    int i;
    // Initialise bgcolor, colors dependent on bgcolor, display strings of numerical vars, tracer directions to positive x
    col_r[101]=28; col_g[101]=33; col_b[101]=47; 
    updateColTheme();
    tweakDispInit();
    fourierDispInit();
    for(i=0;i<100;i++) dir[i]=1;
    // Initialise display strings and filenames for tracer types  
    snprintf (G_frStr, 10, "%lg", G_fr); snprintf (G_ampStr, 10, "%lg", G_amp); snprintf(dxStr, 18, "%.2lf", dx);
    strcpy( imgNames[2], "tracer_3.bmp"); strcpy( imgNames[3], "tracer_4.bmp"); strcpy( imgNames[4], "tracer_5.bmp");
    strcpy( imgNames[5], "tracer_6.bmp"); strcpy( imgNames[6], "tracer_7.bmp"); strcpy( imgNames[7], "tracer_8.bmp");
    strcpy( imgNames[8], "tracer_9.bmp"); strcpy( imgNames[0], "Ball (default)"); strcpy( imgNames[1], "Square bullet");
    strcpy( imgNames[9], "Tangent line"); strcpy( imgNames[10], "    Pacman");

    strcpy( type_GStr[0], "(None)"); strcpy( type_GStr[1], "sin"); strcpy( type_GStr[2], "cos"); strcpy( type_GStr[3], "tan"); 
    strcpy( type_GStr[4], "cosec"); strcpy( type_GStr[5], "sec"); strcpy( type_GStr[6], "cot"); strcpy( type_GStr[7], "arcsin"); 
    strcpy( type_GStr[8], "arccos"); strcpy( type_GStr[9], "arctan"); strcpy( type_GStr[10], "arccosec"); strcpy( type_GStr[11], "arcsec"); 
    strcpy( type_GStr[12], "arccot"); strcpy( type_GStr[13], "sinh"); strcpy( type_GStr[14], "cosh"); strcpy( type_GStr[15], "tanh"); 
    strcpy( type_GStr[16], "cosech"); strcpy( type_GStr[17], "sech"); 
    strcpy( type_GStr[18], "coth"); strcpy( type_GStr[19], "  log_b"); strcpy( type_GStr[20], "b^x"); 

    // Set timers for all animations
    iSetTimer(25, trace);
    iSetTimer(125, pathAnim);
    iSetTimer(200, blinkAnim);

    iInitialize(1360, 700, "Curves (1905084)");
    
    return 0;

};