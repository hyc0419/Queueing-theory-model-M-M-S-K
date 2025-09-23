#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>
#include <limits.h> 
/* Define the constants. */
#define MODLUS 2147483647
#define MULT1 24112
#define MULT2 26143
#define S 2
#define K 8

//隨機數---------------------------------------------------------------------------------------------
/* Set the default seeds for all 100 streams. */
static long zrng[] =
{ 	1,
	1973272912, 281629770, 20006270,1280689831,2096730329,1933576050,
	913566091, 246780520,1363774876, 604901985,1511192140,1259851944,
	824064364, 150493284, 242708531, 75253171,1964472944,1202299975,
	233217322,1911216000, 726370533, 403498145, 993232223,1103205531,
	762430696,1922803170,1385516923, 76271663, 413682397, 726466604,
	336157058,1432650381,1120463904, 595778810, 877722890,1046574445,
	68911991,2088367019, 748545416, 622401386,2122378830, 640690903,
	1774806513,2132545692,2079249579, 78130110, 852776735,1187867272,
	1351423507,1645973084,1997049139, 922510944,2045512870, 898585771,
	243649545,1004818771, 773686062, 403188473, 372279877,1901633463,
	498067494,2087759558, 493157915, 597104727,1530940798,1814496276,
	536444882,1663153658, 855503735, 67784357,1432404475, 619691088,
	119025595, 880802310, 176192644,1116780070, 277854671,1366580350,
	1142483975,2026948561,1053920743, 786262391,1792203830,1494667770,
	1923011392,1433700034,1244184613,1147297105, 539712780,1545929719,
	190641742,1645390429, 264907697, 620389253,1502074852, 927711160,
	364849192,2049576050, 638580085, 547070247 
};

/* Generate the next random number. */
double lcgrand(int stream)
{
	long zi, lowprd, hi31;
	
	zi = zrng[stream];
	lowprd = (zi & 65535) * MULT1;
	hi31 = (zi >> 16) * MULT1 + (lowprd >> 16);
	zi = ((lowprd & 65535) - MODLUS) +
	((hi31 & 32767) << 16) + (hi31 >> 15);
	if (zi < 0) zi += MODLUS;
	lowprd = (zi & 65535) * MULT2;
	hi31 = (zi >> 16) * MULT2 + (lowprd >> 16);
	zi = ((lowprd & 65535) - MODLUS) +
	((hi31 & 32767) << 16) + (hi31 >> 15);
	if (zi < 0) zi += MODLUS;
	zrng[stream] = zi;
	return (zi >> 7 | 1) / 16777216.0;
}

/* Set the current zrng for stream "stream" to zset. */
void lcgrandst (long zset, int stream) 
{
	zrng[stream] = zset;
}

/* Return the current zrng for stream "stream". */
long lcgrandgt (int stream) 
{
	return zrng[stream];
}
//隨機數------------------------------------------------------------------------------------------

struct node{

    double arrivaltime;
    double servicetime;
    double departime;
    struct node *l;

};

struct hd{

    int serv_used;
    double departime;

};

struct node *creat(double arrivaltime , double servicetime){

    struct node *temp=(struct node *)malloc(sizeof(struct node));
    temp->arrivaltime= arrivaltime;
    temp->servicetime= servicetime;
    temp->departime= arrivaltime; //如果沒有被服務的封包會當下就離開
    temp->l=NULL;
    return temp;

}
//封包抵達時間
double poisson_intertime(double rate){

    double interarrivaltime= -log(lcgrand(1))* (60.0/rate);
//    double interarrivaltime= -log( (double)rand() / (double)( RAND_MAX ))* (60.0/rate);
    return interarrivaltime;

}
//檢查目前伺服器最小離開時間
double find_min_depart(struct hd *Serv, int s){

    double min= Serv[0].departime;
    int i;
    for(i=1;i<s;i++){

        if(min > Serv[i].departime){

            min=Serv[i].departime;

        }

    }
    return min;

}
//生成封包資訊
struct node *getarrivalpacketinfo(double lambda, double mu, int n){

    double tempservicetime = poisson_intertime(mu);
    struct node *root = creat(0, tempservicetime);
    struct node *temp=root;
    int i;
    for(i=0;i<n-1;i++){

        double tempintertime = poisson_intertime(lambda);
        tempservicetime = poisson_intertime(mu);
        temp->l= creat(temp->arrivaltime + tempintertime, tempservicetime);
        temp=temp->l;

    }
    return root;
}

int num_Q=0;
struct Queue{

    struct node *l;

};

void arrive(struct hd Serv[], struct node *temp, int s , int k, double mu, struct Queue Q[]){

    int j;
    for(j=0;j<s;j++){

        if(Serv[j].serv_used == 0){
            
            temp->departime=temp->arrivaltime + temp->servicetime;
            Serv[j].departime = temp->departime;
            Serv[j].serv_used = 1;
            return;

        }

    }

    if(num_Q <= k-1){

        Q[num_Q].l=temp;
        num_Q++;

    }

    return;

}

void departure(struct hd Serv[], int s , int k, double mu, struct Queue Q[]){

    int i;
    for(i=0;i<s;i++){

        if(find_min_depart(Serv, s) == Serv[i].departime){

            break;

        }

    }
    if(num_Q == 0){

        Serv[i].serv_used = 0;
        Serv[i].departime = INT_MAX;

    }

    else{

        num_Q--;
        Q->l->departime = Serv[i].departime + Q->l->servicetime;
        Serv[i].departime = Q->l->departime;

        for(i=0;i<num_Q+1;i++){
            
            Q[i].l=Q[i+1].l;
            Q[i+1].l=NULL;

        }

    }

}
//數理分析------------------------------------------------------------------------------------
int factorial(int n){
    int i;
    if(n==1 || n==0){
        return 1;
    }
    for(i=n-1;i>0;i--){
        n=n*i;
    }
    return n;
}

double P0(float lambda, float mu, int s, int k){
	double sum = 1.0;
	int i;
	for(i=1;i<s;i++){
        sum += pow((lambda/mu),i) / factorial(i);
    }
	sum += (pow((lambda/mu), s) / factorial(s) * ((1-pow((lambda/(s*mu)), (k-s+1))) / (1-(lambda/(s*mu)))));
	return 1/sum;
}

double Pn(float lambda, float mu, int s, int k, int n){
	if(n<s){
        return ((pow((lambda/mu), n)) / factorial(n)) * P0(lambda, mu, s, k);
    }
		
	return ((pow((lambda/mu), n)) / (factorial(s) * pow(s, (n-s)))) * P0(lambda, mu, s, k);
}

double Lq_math(float lambda, float mu, int s, int k){
	double result = 0;
	int i;
	for(i=s;i<k+1;i++)
		result += (i-s)*Pn(lambda, mu, s, k, i);
	return result;
}

double L_math(float lambda, float mu, int s, int k){
	double sum1 = 0, sum2 = 0, temp;
	int i;
	for(i=0;i<s;i++){

		temp = Pn(lambda, mu, s, k, i);
		sum1 += i*temp;
		sum2 += temp; 

	}
	sum2 = s * (1- sum2);
	return (Lq_math(lambda, mu, s, k) + sum1 + sum2);
}
//數理分析------------------------------------------------------------------------------------


int main(){

    FILE *fp;
    fp = fopen("mmsk.csv", "a+"); 
    double W = 0.0, Wq =0.0, L = 0.0, Lq= 0.0;
	int n=1000000, i,j;
    double lambda=25.0;
    double mu=40.0;

        struct hd Serv[S];
        for(i=0;i<S;i++){

            Serv[i].serv_used=0;
            Serv[i].departime=INT_MAX; //get huge value
        
        }

        struct Queue Q[K+1];
        for(i=0;i<K;i++){

            Q[i].l=NULL;

        }
		//產生n個封包的抵達時間和服務時間
        struct node *root=getarrivalpacketinfo(lambda, mu, n);
        struct node *temp=root;

		//當temp不為NULL時，進行封包抵達和離開的處理
        while(temp != NULL){
        
            if(temp->arrivaltime < find_min_depart(Serv, S)){

                arrive(Serv, temp, S , K, mu, Q);
                temp=temp->l;

            }

            else{

                departure(Serv, S, K, mu, Q);

            }
        
        }

        temp=root;
        int count=0;
        double totalservicetime=0, sum=0.0;

		//計算L, Lq, W, Wq
        while(temp->l!= NULL){

            if(temp->departime != temp->arrivaltime){
                
                W += temp->departime - temp->arrivaltime;
                totalservicetime += temp->servicetime;
                count++;

            }
            temp=temp->l;

        }


        printf("%.8lf, ", W/(temp->departime)); //L
        printf("%.8lf, ", (W-totalservicetime)/(temp->departime)); //Lq
        printf("%.8lf, ", W/count/60); //W
        printf("%.8lf\n", (W-totalservicetime)/count/60); //Wq

        printf("%.8lf, ", L_math(lambda, mu, S, K)); //L
	    printf("%.8lf, ", Lq_math(lambda, mu, S, K)); //Lq
	    printf("%.8lf, ", L_math(lambda, mu, S, K) / (lambda*(1-Pn(lambda, mu, S, K, K)))); //W
	    printf("%.8lf ", Lq_math(lambda, mu, S, K) / (lambda*(1-Pn(lambda, mu, S, K, K)))); //Wq



        fprintf(fp,"%.8lf,", W/(temp->departime)); //L
        fprintf(fp,"%.8lf,", (W-totalservicetime)/(temp->departime));//Lq
        fprintf(fp,"%.8lf,", W/count/60);//W
        fprintf(fp,"%.8lf,,", (W-totalservicetime)/count/60);//Wq
    
        fprintf(fp,"%.8lf,", L_math(lambda, mu, S, K));//L
	    fprintf(fp,"%.8lf,", Lq_math(lambda, mu, S, K));//Lq
	    fprintf(fp,"%.8lf,", L_math(lambda, mu, S, K) / (lambda*(1-Pn(lambda, mu, S, K, K))));//W
	    fprintf(fp,"%.8lf\n", Lq_math(lambda, mu, S, K) / (lambda*(1-Pn(lambda, mu, S, K, K))));//Wq

    fclose(fp);
    return 0;
}
