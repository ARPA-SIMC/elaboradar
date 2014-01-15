#define  ANAP_OK  0
#define  ANAP_YES 1
#define  ANAP_NODAT  2
#define  ANAP_NOCONTROL  3
#define  BBMAX 50
#define  DTLIM 4.0
#define  DTMAX 6.0
#define  DRLIM 100000 //cambiato
#define  ERRPT 0.1
#define  BETA  0.00732 //(parametro climatico adj calcolato da ROBBY)

#ifdef __cplusplus
extern "C" {
#endif
float qBB(unsigned char bbc,float dr,float dt);
float qCl( unsigned char  clc);
float qDist( float dr);
float qVol( float dh, float dhst);
float qAtt(float PIA);
float qVpr(float  dZ,float sdevZ);
#ifdef __cplusplus
}
#endif
