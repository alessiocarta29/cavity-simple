//CAVITA' QUADRATA LID-DRIVEN CON VELOCITA' COSTANTE U
//METODO FV CON ALGORITMO SIMPLE, GRIGLIA CARTESIANA UNIFORME COLLOCATA,
//DISCRETIZZAZIONE CDS, SISTEMI LINEARI CON SIP, PORTATE CON RHIE-CHOW

#include<cstdio>
#include<cmath>
#include<vector>
#include<algorithm>
#include<ctime>

//TUTTI I CAMPI SONO VETTORI DI M = N*N ELEMENTI, CON k = j + N*i
//i = RIGA (DIREZIONE y), j = COLONNA (DIREZIONE x)
//E = k+1,  W = k-1,  N = k+N,  S = k-N


//==========================  STRUTTURE  ==========================

//PARAMETRI DELLA SIMULAZIONE, RACCOLTI IN UN UNICO POSTO COSI' LE
//FUNZIONI NE RICEVONO UNO SOLO INVECE DI DIECI ARGOMENTI SCIOLTI
struct Parametri
{
    int N;            //CELLE PER DIREZIONE
    double L;         //LATO DELLA CAVITA'
    double U;         //VELOCITA' DEL COPERCHIO
    double Re;
    double rho;
    double au;        //SOTTORILASSAMENTO QUANTITA' DI MOTO
    double ap;        //SOTTORILASSAMENTO PRESSIONE
    double gamma;     //1 = CDS, 0 = UDS, IN MEZZO UNA MISCELA DEI DUE
    int nswUV;        //SPAZZATE SIP PER LA QUANTITA' DI MOTO
    int nswP;         //SPAZZATE SIP PER LA PRESSIONE
    double tol;       //RIDUZIONE DEL RESIDUO RICHIESTA (NON UN VALORE ASSOLUTO)
    int maxIter;
    int verbose;      //1 = STAMPA UNA RIGA PER ITERAZIONE, 0 = ZITTO
};

//I CINQUE COEFFICENTI E IL TERMINE NOTO DI UN SISTEMA A 5 DIAGONALI
struct Sistema
{
    std::vector<double> aE, aW, aN, aS, aP, q;
};

//COSA RESTITUISCE UNA SIMULAZIONE
struct Risultato
{
    int iter;
    double erru, errv, errm, err;
    int convergiuto;   //1 SE HA RAGGIUNTO tol, 0 SE E' USCITO PER maxIter
    double tempo;      //SECONDI
};


//==========================  UTILITA'  ==========================

//VALORI DI PARTENZA
Parametri parametriDefault()
{
    Parametri par;
    par.N = 128;
    par.L = 1.0;
    par.U = 5.0;
    par.Re = 1000.0;
    par.rho = 1000.0;
    par.au = 0.8;
    par.ap = 0.2;
    par.gamma = 1.0;      //CDS PURO
    par.nswUV = 3;
    par.nswP = 20;
    //tol E' UNA RIDUZIONE RELATIVA: 1e-4 VUOL DIRE QUATTRO ORDINI DI
    //GRANDEZZA IN MENO RISPETTO AL RESIDUO DI PARTENZA. FERZIGER 5.7
    //CONSIGLIA DA TRE A CINQUE ORDINI PER LE ITERAZIONI ESTERNE
    par.tol = 1e-4;
    par.maxIter = 5000;
    par.verbose = 1;
    return par;
}

//GRANDEZZE RICAVATE DAGLI ALTRI PARAMETRI: NON SI MEMORIZZANO, SI
//RICALCOLANO QUANDO SERVONO, COSI' NON POSSONO RESTARE DISALLINEATE
double passo(Parametri par)
{
    return par.L/par.N;
}

double viscosita(Parametri par)
{
    return par.rho*par.U*par.L/par.Re;
}

//ALLOCA I SEI VETTORI DI UN SISTEMA E LI AZZERA.
//assign(M,0.0) PORTA IL VETTORE A M ELEMENTI TUTTI UGUALI A ZERO
void dimensiona(Sistema &s, int M)
{
    s.aE.assign(M,0.0);
    s.aW.assign(M,0.0);
    s.aN.assign(M,0.0);
    s.aS.assign(M,0.0);
    s.aP.assign(M,0.0);
    s.q.assign(M,0.0);
}


//====================  SOLUTORE LINEARE (SIP)  ====================

//METODO SIP DI STONE PER MATRICI A 5 DIAGONALI
//RISOLVE:  aP[k]*f[k] + aE[k]*f[k+1] + aW[k]*f[k-1] + aN[k]*f[k+N] + aS[k]*f[k-N] = q[k]
//nsw = NUMERO DI SPAZZATE. RESTITUISCE IL RESIDUO INIZIALE.
//FUNZIONA SOLO SE aP>0 E I COEFFICENTI DEI VICINI SONO <=0
double SIP(int N, Sistema &s, std::vector<double> &f, int nsw)
{
    int M = N*N;
    double alfa = 0.92;

    std::vector<double> LW(M,0.0), LS(M,0.0), LP(M,0.0);
    std::vector<double> UE(M,0.0), UN(M,0.0), R(M,0.0);

    //FATTORIZZAZIONE INCOMPLETA A = L*U
    //L HA LE DIAGONALI LS(k-N), LW(k-1), LP(k)
    //U HA LE DIAGONALI 1(k), UE(k+1), UN(k+N)
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;

            double UEs = 0.0, UNs = 0.0, UEw = 0.0, UNw = 0.0;
            if(i!=0) { UEs = UE[k-N]; UNs = UN[k-N]; }
            if(j!=0) { UEw = UE[k-1]; UNw = UN[k-1]; }

            LS[k] = s.aS[k]/(1.0 + alfa*UEs);
            LW[k] = s.aW[k]/(1.0 + alfa*UNw);

            double p1 = alfa*LS[k]*UEs;
            double p2 = alfa*LW[k]*UNw;

            LP[k] = 1.0/(s.aP[k] + p1 + p2 - LS[k]*UNs - LW[k]*UEw);   //LP = 1/L^P

            UE[k] = (s.aE[k] - p1)*LP[k];
            UN[k] = (s.aN[k] - p2)*LP[k];
        }
    }

    double res0 = 0.0;

    for(int it=0;it<nsw;it++)
    {
        //RESIDUO r = q - A*f  E DISCESA IN AVANTI  R = L^-1 * r
        double res = 0.0;
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;

                double r = s.q[k] - s.aP[k]*f[k];
                if(j!=N-1) r = r - s.aE[k]*f[k+1];
                if(j!=0)   r = r - s.aW[k]*f[k-1];
                if(i!=N-1) r = r - s.aN[k]*f[k+N];
                if(i!=0)   r = r - s.aS[k]*f[k-N];

                res = res + std::abs(r);

                if(i!=0) r = r - LS[k]*R[k-N];
                if(j!=0) r = r - LW[k]*R[k-1];

                R[k] = r*LP[k];
            }
        }
        if(it==0) res0 = res;

        //RISALITA ALL'INDIETRO  R = U^-1 * R  E CORREZIONE DELLA VARIABILE
        for(int i=N-1;i>=0;i--)
        {
            for(int j=N-1;j>=0;j--)
            {
                int k = j + N*i;

                if(i!=N-1) R[k] = R[k] - UN[k]*R[k+N];
                if(j!=N-1) R[k] = R[k] - UE[k]*R[k+1];

                f[k] = f[k] + R[k];
            }
        }
    }

    return res0;
}


//========================  ASSEMBLAGGIO  ========================

//VALORE DI FACCIA CON LO SCHEMA UPWIND: SI PRENDE LA CELLA DA CUI IL
//FLUIDO ARRIVA. m E' LA PORTATA CONTATA POSITIVA DA "prima" VERSO "dopo"
double upwind(double m, double fiPrima, double fiDopo)
{
    if(m>0.0) { return fiPrima; }
    return fiDopo;
}

//COEFFICENTI DI UNA DELLE DUE EQUAZIONI DI QUANTITA' DI MOTO.
//direzione = 0 -> COMPONENTE x,  direzione = 1 -> COMPONENTE y
//fiParete E' LA VELOCITA' DEL COPERCHIO PER QUESTA COMPONENTE:
//U PER LA x, 0 PER LA y, PERCHE' IL COPERCHIO SCORRE SOLO IN ORIZZONTALE
void coefficientiQdm(Parametri par, int direzione, double fiParete,
                     std::vector<double> &fi, std::vector<double> &p,
                     std::vector<double> &fe, std::vector<double> &fn,
                     Sistema &s)
{
    int N = par.N;
    double dx = passo(par);
    double mu = viscosita(par);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;

            //LA PORTATA DI FACCIA NON SI INTERPOLA QUI: SI LEGGE DA fe/fn,
            //COSTRUITE CON RHIE-CHOW DOPO LA QUANTITA' DI MOTO.
            //LA FACCIA OVEST DELLA CELLA k E' LA FACCIA EST DELLA CELLA k-1
            //
            //NELLA MATRICE VA SEMPRE L'UPWIND, MAI IL CDS: SOLO L'UPWIND
            //SODDISFA aP >= somma|a_nb|, CHE E' LA CONDIZIONE SUFFICIENTE
            //PERCHE' I SOLUTORI ITERATIVI CONVERGANO (FERZIGER 5.8).
            //min(m,0) E' LA PORTATA SOLO SE ENTRA DAL VICINO, max(m,0)
            //SOLO SE ESCE VERSO DI LUI: E' L'UPWIND SCRITTO SENZA if
            if(j!=N-1) { s.aE[k] =  std::min(fe[k],0.0)    - mu; }
            if(j!=0)   { s.aW[k] = -std::max(fe[k-1],0.0)  - mu; }
            if(i!=N-1) { s.aN[k] =  std::min(fn[k],0.0)    - mu; }
            if(i!=0)   { s.aS[k] = -std::max(fn[k-N],0.0)  - mu; }

            //PRESSIONE SULLE FACCE. SE LA CELLA E' A PARETE SI USA
            //LA PRESSIONE DELLA CELLA STESSA
            double pW = p[k], pE = p[k], pS = p[k], pN = p[k];
            if(j!=0)   { pW = p[k-1]; }
            if(j!=N-1) { pE = p[k+1]; }
            if(i!=0)   { pS = p[k-N]; }
            if(i!=N-1) { pN = p[k+N]; }

            if(direzione==0) { s.q[k] = dx*0.5*(pW - pE); }
            else             { s.q[k] = dx*0.5*(pS - pN); }

            //NO-SLIP. LA PARETE STA A dx/2 DAL CENTRO CELLA, NON A dx,
            //QUINDI IL SUO COEFFICENTE DIFFUSIVO E' mu*dx/(dx/2) = 2*mu.
            //ATTRAVERSO LA PARETE NON PASSA MASSA, QUINDI NIENTE CONVEZIONE.
            //LA PARTE INCOGNITA DEL FLUSSO VA SU aP, QUELLA NOTA SU q.
            double aWall = 0.0;
            if(j==N-1) { aWall = aWall + 2*mu; }
            if(j==0)   { aWall = aWall + 2*mu; }
            if(i==0)   { aWall = aWall + 2*mu; }
            if(i==N-1) { aWall = aWall + 2*mu; s.q[k] = s.q[k] + 2*mu*fiParete; }

            //PORTATA NETTA USCENTE DALLA CELLA, SERVE PER aP
            double dm = 0.0;
            if(j!=N-1) { dm = dm + fe[k];   }
            if(j!=0)   { dm = dm - fe[k-1]; }
            if(i!=N-1) { dm = dm + fn[k];   }
            if(i!=0)   { dm = dm - fn[k-N]; }

            //aP E' MENO LA SOMMA DEI VICINI PIU' LA PORTATA NETTA.
            //I VICINI CHE NON ESISTONO VALGONO 0 E NON DANNO CONTRIBUTO
            s.aP[k] = -(s.aE[k] + s.aW[k] + s.aN[k] + s.aS[k]) + dm + aWall;

            //CORREZZIONE DIFFERITA. LA MATRICE E' UPWIND, MA LO SCHEMA CHE
            //SI VUOLE E' gamma*CDS + (1-gamma)*UDS. LA DIFFERENZA TRA I DUE
            //SI METTE NEL TERMINE NOTO CALCOLATA COL CAMPO DELL'ITERAZIONE
            //PRECEDENTE. A CONVERGENZA fi NON CAMBIA PIU', L'UPWIND IMPLICITO
            //SI CANCELLA CON QUELLO ESPLICITO E RESTA LO SCHEMA VOLUTO,
            //MENTRE LA MATRICE RESTA QUELLA STABILE DELL'UPWIND.
            //dc E' LA SOMMA SULLE QUATTRO FACCE DI m*(fi_UDS - fi_CDS),
            //COL SEGNO POSITIVO SE LA PORTATA ESCE DALLA CELLA
            double dc = 0.0;
            if(j!=N-1) { dc = dc + fe[k]  *( upwind(fe[k],  fi[k],  fi[k+1]) - 0.5*(fi[k]+fi[k+1]) ); }
            if(j!=0)   { dc = dc - fe[k-1]*( upwind(fe[k-1],fi[k-1],fi[k])   - 0.5*(fi[k-1]+fi[k]) ); }
            if(i!=N-1) { dc = dc + fn[k]  *( upwind(fn[k],  fi[k],  fi[k+N]) - 0.5*(fi[k]+fi[k+N]) ); }
            if(i!=0)   { dc = dc - fn[k-N]*( upwind(fn[k-N],fi[k-N],fi[k])   - 0.5*(fi[k-N]+fi[k]) ); }
            s.q[k] = s.q[k] + par.gamma*dc;

            //SOTTORILASSAMENTO IMPLICITO
            s.aP[k] = s.aP[k]/par.au;
            s.q[k] = s.q[k] + (1-par.au)*s.aP[k]*fi[k];
        }
    }
}

//PORTATE SULLE FACCE CON INTERPOLAZIONE DI QUANTITA' DI MOTO (RHIE-CHOW).
//LA MEDIA (u_P+u_E)/2 SI PORTA DIETRO I GRADIENTI CENTRATI DELLE DUE CELLE,
//NEI QUALI (p_E - p_P) SI SEMPLIFICA: LA PORTATA NON VEDREBBE IL SALTO DI
//PRESSIONE ATTRAVERSO LA FACCIA CHE STA ATTRAVERSANDO. LA PARENTESI TOGLIE
//IL GRADIENTE MEDIATO E RIMETTE QUELLO A DUE PUNTI. SU CAMPO REGOLARE I DUE
//QUASI COINCIDONO, SU UNA SCACCHIERA NO E LA SMORZANO. COSI' LA PORTATA
//RISPONDE A (p_E - p_P) CON COEFFICENTE (V/aP)_e, CHE E' ESATTAMENTE QUELLO
//SCRITTO IN aEp: LE DUE EQUAZIONI SONO COERENTI
void portateFacce(Parametri par,
                  std::vector<double> &ustar, std::vector<double> &vstar,
                  std::vector<double> &p,
                  std::vector<double> &aPu, std::vector<double> &aPv,
                  std::vector<double> &fe, std::vector<double> &fn)
{
    int N = par.N;
    int M = N*N;
    double dx = passo(par);

    //GRADIENTI DI p AL CENTRO CELLA (STESSE FACCE USATE PER q)
    std::vector<double> gpx(M,0.0), gpy(M,0.0);
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            double pW = p[k], pE = p[k], pS = p[k], pN = p[k];
            if(j!=0)   { pW = p[k-1]; }
            if(j!=N-1) { pE = p[k+1]; }
            if(i!=0)   { pS = p[k-N]; }
            if(i!=N-1) { pN = p[k+N]; }
            gpx[k] = (pE - pW)*0.5/dx;
            gpy[k] = (pN - pS)*0.5/dx;
        }
    }

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            if(j!=N-1)
            {
            double C = 0.5*dx*dx*(1.0/aPu[k] + 1.0/aPu[k+1]);
            fe[k] = par.rho*dx*( 0.5*(ustar[k] + ustar[k+1])
                  - C*((p[k+1] - p[k])/dx - 0.5*(gpx[k] + gpx[k+1])) );
            }
            if(i!=N-1)
            {
            double C = 0.5*dx*dx*(1.0/aPv[k] + 1.0/aPv[k+N]);
            fn[k] = par.rho*dx*( 0.5*(vstar[k] + vstar[k+N])
                  - C*((p[k+N] - p[k])/dx - 0.5*(gpy[k] + gpy[k+N])) );
            }
        }
    }
}

//PORTATA MASSICA IN ECCEDENZA DA OGNI CELLA
//LE FACCE DI PARETE NON COMPAIONO PERCHE' LI' LA PORTATA E' NULLA
void sbilancio(Parametri par, std::vector<double> &fe, std::vector<double> &fn,
               std::vector<double> &dmp)
{
    int N = par.N;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            double me = 0.0, mw = 0.0, mn = 0.0, ms = 0.0;
            if(j!=N-1) { me = fe[k];   }
            if(j!=0)   { mw = fe[k-1]; }
            if(i!=N-1) { mn = fn[k];   }
            if(i!=0)   { ms = fn[k-N]; }
            dmp[k] = mn - ms + me - mw;
        }
    }
}

//COEFFICENTI DELL'EQUAZIONE DELLA CORREZZIONE DI PRESSIONE.
//IL COEFFICENTE E' -rho*(V/aP) DELLA FACCIA, E LA FACCIA E' UNA SOLA:
//SI USA LA MEDIA DELLE DUE CELLE, COSI' aE[k] E aW[k+1] COINCIDONO.
//ATTRAVERSO UNA PARETE NON C'E' CORREZZIONE, QUINDI IL COEFFICENTE E' 0
void coefficientiPressione(Parametri par,
                           std::vector<double> &aPu, std::vector<double> &aPv,
                           std::vector<double> &dmp, Sistema &s)
{
    int N = par.N;
    int M = N*N;
    double dx = passo(par);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            if(j!=N-1) { s.aE[k] = -par.rho*0.5*dx*dx*(1.0/aPu[k] + 1.0/aPu[k+1]); }
            if(j!=0)   { s.aW[k] = -par.rho*0.5*dx*dx*(1.0/aPu[k] + 1.0/aPu[k-1]); }
            if(i!=N-1) { s.aN[k] = -par.rho*0.5*dx*dx*(1.0/aPv[k] + 1.0/aPv[k+N]); }
            if(i!=0)   { s.aS[k] = -par.rho*0.5*dx*dx*(1.0/aPv[k] + 1.0/aPv[k-N]); }
            s.aP[k] = -(s.aE[k] + s.aW[k] + s.aN[k] + s.aS[k]);
        }
    }

    //IL TERMINE NOTO E' MENO LA PORTATA IN ECCEDENZA
    for(int k=0;k<M;k++)
    {
        s.q[k] = -dmp[k];
    }

    //LA PRESSIONE E' DEFINITA A MENO DI UNA COSTANTE, QUINDI LA MATRICE E'
    //SINGOLARE: SI FISSA p' = 0 IN UNA CELLA DI RIFERIMENTO (LA 0)
    s.aE[0] = 0.0;
    s.aW[0] = 0.0;
    s.aN[0] = 0.0;
    s.aS[0] = 0.0;
    s.aP[0] = 1.0;
    s.q[0] = 0.0;
}


//=========================  CORREZZIONI  =========================

//CORREZZIONE DELLE PORTATE DI FACCIA.
//m' = -rho*(V/aP)*(p'_E - p'_P), CHE E' GIA' IL COEFFICENTE aE DEL SISTEMA
//DELLA PRESSIONE: NON SERVE RICALCOLARLO. E' IL PASSO CHE RENDE dmp NULLO
void correggiPortate(Parametri par, std::vector<double> &pc, Sistema &sp,
                     std::vector<double> &fe, std::vector<double> &fn)
{
    int N = par.N;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            if(j!=N-1) { fe[k] = fe[k] + sp.aE[k]*(pc[k+1] - pc[k]); }
            if(i!=N-1) { fn[k] = fn[k] + sp.aN[k]*(pc[k+N] - pc[k]); }
        }
    }
}

//CAMPO DELLE VELOCITA' CORRETTE E PRESSIONE SOTTORILASSATA.
//STESSA FORMULA DELLE FACCE MA COL GRADIENTE CENTRATO SULLA CELLA.
//SULLA PARETE p' DI FACCIA E' QUELLO DELLA CELLA, COME PER p
void correggiCelle(Parametri par, std::vector<double> &pc,
                   std::vector<double> &aPu, std::vector<double> &aPv,
                   std::vector<double> &ustar, std::vector<double> &vstar,
                   std::vector<double> &u, std::vector<double> &v,
                   std::vector<double> &p)
{
    int N = par.N;
    double dx = passo(par);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            double pcW = pc[k], pcE = pc[k], pcS = pc[k], pcN = pc[k];
            if(j!=0)   { pcW = pc[k-1]; }
            if(j!=N-1) { pcE = pc[k+1]; }
            if(i!=0)   { pcS = pc[k-N]; }
            if(i!=N-1) { pcN = pc[k+N]; }

            //LA VELOCITA' SI CORREGGE IN PIENO, LA PRESSIONE SOTTORILASSATA
            u[k] = ustar[k] + dx*0.5*(pcW - pcE)/aPu[k];
            v[k] = vstar[k] + dx*0.5*(pcS - pcN)/aPv[k];
            p[k] = p[k] + pc[k]*par.ap;
        }
    }
}


//========================  CICLO SIMPLE  ========================

//RISOLVE IL PROBLEMA CON I PARAMETRI DATI. u, v, p ENTRANO COME STIMA
//INIZIALE ED ESCONO CON LA SOLUZIONE
Risultato risolviCavita(Parametri par, std::vector<double> &u,
                        std::vector<double> &v, std::vector<double> &p)
{
    int N = par.N;
    int M = N*N;

    Sistema su, sv, sp;      //QUANTITA' DI MOTO x, y E PRESSIONE
    dimensiona(su,M);
    dimensiona(sv,M);
    dimensiona(sp,M);

    //PORTATE MEMORIZZATE SULLA FACCIA EST E NORD DI OGNI CELLA.
    //NON SI RICALCOLANO DA ZERO OGNI VOLTA: SONO UNO STATO CHE SOPRAVVIVE
    //TRA UN'ITERAZIONE E L'ALTRA, ALTRIMENTI LA CORREZZIONE DI p' VA PERSA
    std::vector<double> fe(M,0.0), fn(M,0.0), dmp(M,0.0);
    std::vector<double> ustar(M,0.0), vstar(M,0.0), pc(M,0.0);

    Risultato ris;
    ris.iter = 0;
    ris.err = 0.0;
    ris.convergiuto = 0;

    //RESIDUI DI RIFERIMENTO, UNO PER EQUAZIONE. SI PRENDONO ALLA PRIMA
    //ITERAZIONE E POI SI CONFRONTA IL RAPPORTO res/res0 CON tol.
    //SENZA QUESTO LO STESSO tol IMPORREBBE UNA CONDIZIONE DIVERSA A OGNI N,
    //E I CONFRONTI DI NUMERO DI ITERAZIONI TRA CONFIGURAZIONI DIVERSE
    //MISUREREBBERO IL CRITERIO INVECE DEL SOLUTORE.
    //ALLA PRIMA ITERAZIONE res_v E' ESATTAMENTE ZERO (v E' NULLA OVUNQUE E
    //IL COPERCHIO NON SPINGE IN VERTICALE), QUINDI SI PRENDE COME
    //RIFERIMENTO IL PRIMO VALORE NON NULLO
    double res0u = 0.0, res0v = 0.0, res0m = 0.0;

    //clock() CONTA IL TEMPO DI CALCOLO; DIVISO PER CLOCKS_PER_SEC DA' I SECONDI
    clock_t t0 = clock();

    do
    {
        coefficientiQdm(par,0,par.U,u,p,fe,fn,su);
        coefficientiQdm(par,1,0.0,  v,p,fe,fn,sv);

        //SI PARTE DAL CAMPO DELL'ITERAZIONE PRECEDENTE COME STIMA INIZIALE.
        //POCHE SPAZZATE: DENTRO AL SIMPLE NON SERVE RISOLVERE A CONVERGENZA
        for(int k=0;k<M;k++) { ustar[k] = u[k]; }
        for(int k=0;k<M;k++) { vstar[k] = v[k]; }
        double resu = SIP(N,su,ustar,par.nswUV);
        double resv = SIP(N,sv,vstar,par.nswUV);

        portateFacce(par,ustar,vstar,p,su.aP,sv.aP,fe,fn);
        sbilancio(par,fe,fn,dmp);

        coefficientiPressione(par,su.aP,sv.aP,dmp,sp);

        //SI RIPARTE SEMPRE DA p' = 0.
        //QUI SERVONO PIU' SPAZZATE: E' L'EQUAZIONE ELLITTICA
        for(int k=0;k<M;k++) { pc[k] = 0.0; }
        SIP(N,sp,pc,par.nswP);

        correggiPortate(par,pc,sp,fe,fn);
        correggiCelle(par,pc,su.aP,sv.aP,ustar,vstar,u,v,p);

        //CRITERIO DI ARRESTO SUI RESIDUI.
        //resu E resv LI RESTITUISCE GIA' IL SIP: SONO LA SOMMA DI |q - A*f|
        //CALCOLATA PRIMA DI RISOLVERE, CIOE' DI QUANTO L'EQUAZIONE DISCRETA
        //NON E' SODDISFATTA DAL CAMPO CORRENTE.
        //resm E' IL RESIDUO DELLA CONTINUITA': LA SOMMA DEGLI SBILANCI
        double resm = 0.0;
        for(int k=0;k<M;k++)
        {
            resm = resm + std::abs(dmp[k]);
        }

        //IL PRIMO VALORE NON NULLO DIVENTA IL RIFERIMENTO DI QUELL'EQUAZIONE
        if(res0u==0.0) { res0u = resu; }
        if(res0v==0.0) { res0v = resv; }
        if(res0m==0.0) { res0m = resm; }

        //DI QUANTO IL RESIDUO E' SCESO RISPETTO ALLA PARTENZA
        ris.erru = 1.0;
        ris.errv = 1.0;
        ris.errm = 1.0;
        if(res0u>0.0) { ris.erru = resu/res0u; }
        if(res0v>0.0) { ris.errv = resv/res0v; }
        if(res0m>0.0) { ris.errm = resm/res0m; }

        //SI USA IL PIU' GRANDE DEI TRE
        ris.err = std::max(ris.erru,ris.errv);
        ris.err = std::max(ris.err,ris.errm);

        ris.iter = ris.iter + 1;
        if(par.verbose==1)
        {
            printf("iter %5d   res_u %10.3e   res_v %10.3e   res_massa %10.3e\n",
                   ris.iter,ris.erru,ris.errv,ris.errm);
        }

        //TETTO ALLE ITERAZIONI PER NON RESTARE BLOCCATI SE LA SOLUZIONE DIVERGE
    }while(ris.err>par.tol && ris.iter<par.maxIter);

    ris.tempo = (double)(clock() - t0)/CLOCKS_PER_SEC;
    if(ris.err<=par.tol) { ris.convergiuto = 1; }

    return ris;
}


//======================  INTENSITA' DEI VORTICI  ======================

//VORTICITA' AL CENTRO CELLA: omega = dv/dx - du/dy.
//LE DERIVATE SI FANNO CON I VALORI SULLE DUE FACCE OPPOSTE, CHE DISTANO
//ESATTAMENTE dx, QUINDI E' UNA DIFFERENZA CENTRATA. SULLA PARETE NON SI
//INTERPOLA: SI USA IL VALORE VERO, CHE E' NOTO
void vorticita(Parametri par, std::vector<double> &u, std::vector<double> &v,
               std::vector<double> &om)
{
    int N = par.N;
    double dx = passo(par);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;

            double ve = 0.0, vw = 0.0, un = 0.0, us = 0.0;
            if(j!=N-1) { ve = 0.5*(v[k] + v[k+1]); }    //ALTRIMENTI PARETE, v = 0
            if(j!=0)   { vw = 0.5*(v[k] + v[k-1]); }
            if(i!=N-1) { un = 0.5*(u[k] + u[k+N]); }
            else       { un = par.U; }                  //COPERCHIO
            if(i!=0)   { us = 0.5*(u[k] + u[k-N]); }    //ALTRIMENTI FONDO, u = 0

            om[k] = (ve - vw)/dx - (un - us)/dx;
        }
    }
}

//FUNZIONE DI CORRENTE, DEFINITA DA u = dpsi/dy E v = -dpsi/dx.
//SOSTITUENDO NELLA VORTICITA' SI OTTIENE UNA POISSON: lap(psi) = -omega.
//SI RISOLVE CON LO STESSO SIP DEL RESTO DEL PROGRAMMA.
//LE QUATTRO PARETI SONO LINEE DI CORRENTE DI UNA CAVITA' CHIUSA, QUINDI
//psi = 0 SU TUTTE, ED E' UNA DIRICHLET COME QUELLA DELLA VELOCITA':
//PARETE A dx/2, COEFFICENTE 2 INVECE DI 1
//RESTITUISCE DI QUANTO E' SCESO IL RESIDUO: SERVE A SAPERE SE PSI E'
//DAVVERO RISOLTA O SE IL NUMERO CHE SI LEGGE E' SOLO UN'ITERAZIONE A META'
double funzioneCorrente(Parametri par, std::vector<double> &om,
                        std::vector<double> &psi)
{
    int N = par.N;
    int M = N*N;
    double dx = passo(par);

    Sistema s;
    dimensiona(s,M);

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            double conta = 0.0;

            if(j!=N-1) { s.aE[k] = -1.0; conta = conta + 1.0; } else { conta = conta + 2.0; }
            if(j!=0)   { s.aW[k] = -1.0; conta = conta + 1.0; } else { conta = conta + 2.0; }
            if(i!=N-1) { s.aN[k] = -1.0; conta = conta + 1.0; } else { conta = conta + 2.0; }
            if(i!=0)   { s.aS[k] = -1.0; conta = conta + 1.0; } else { conta = conta + 2.0; }

            s.aP[k] = conta;
            s.q[k] = om[k]*dx*dx;
        }
    }

    //QUI SI RISOLVE SUL SERIO, NON DENTRO A UN CICLO ESTERNO.
    //IL NUMERO DI SPAZZATE NON SI PUO' FISSARE A MANO: SU GRIGLIA FINE NE
    //SERVONO MOLTE DI PIU' (A N=128 NE VOGLIONO CIRCA 1500 CONTRO LE 400
    //CHE BASTANO A N=64), E UNA POISSON LASCIATA A META' SPOSTA psi DI UN
    //PAIO DI PUNTI PERCENTUALE SENZA DARE NESSUN SEGNALE DI ERRORE.
    //SI ITERA A BLOCCHI FINCHE' IL RESIDUO NON E' SCESO DI OTTO ORDINI
    for(int k=0;k<M;k++) { psi[k] = 0.0; }

    double primo = 0.0;
    double res = 0.0;
    for(int blocco=0;blocco<300;blocco++)
    {
        res = SIP(N,s,psi,20);
        if(blocco==0) { primo = res; }
        if(res < primo*1e-8) { break; }
    }

    if(primo>0.0) { return res/primo; }
    return 0.0;
}

//CONTROLLO DEL POST-PROCESSING. psi SI PUO' OTTENERE ANCHE INTEGRANDO
//DIRETTAMENTE u LUNGO y (OPPURE -v LUNGO x) PARTENDO DALLA PARETE DOVE
//VALE ZERO. LE DUE STRADE DEVONO DARE LO STESSO RISULTATO: SE NON LO
//DANNO, L'ERRORE E' NEL POST-PROCESSING, NON NEL SOLUTORE.
//verso = 0 -> INTEGRA u LUNGO y,  verso = 1 -> INTEGRA -v LUNGO x
void correnteIntegrata(Parametri par, std::vector<double> &u,
                       std::vector<double> &v, int verso,
                       std::vector<double> &psi)
{
    int N = par.N;
    double dx = passo(par);

    if(verso==0)
    {
        //PER OGNI COLONNA SI RISALE DAL FONDO SOMMANDO u*dx.
        //IL VALORE ALLE FACCE VIENE DALLA SOMMA, QUELLO AL CENTRO CELLA
        //E' LA MEDIA DELLE DUE FACCE CHE LA DELIMITANO
        for(int j=0;j<N;j++)
        {
            double sotto = 0.0;
            for(int i=0;i<N;i++)
            {
                int k = j + N*i;
                double sopra = sotto + u[k]*dx;
                psi[k] = 0.5*(sotto + sopra);
                sotto = sopra;
            }
        }
    }
    else
    {
        for(int i=0;i<N;i++)
        {
            double sinistra = 0.0;
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                double destra = sinistra - v[k]*dx;
                psi[k] = 0.5*(sinistra + destra);
                sinistra = destra;
            }
        }
    }
}

//INTENSITA' DEI VORTICI: SONO GLI ESTREMI DELLA FUNZIONE DI CORRENTE,
//NORMALIZZATI CON U*L. psiMin E' IL VORTICE PRIMARIO, psiMax QUELLI
//SECONDARI NEGLI ANGOLI IN BASSO. psiMax E' MOLTO PIU' SENSIBILE ALLA
//GRIGLIA E ALLO SCHEMA CONVETTIVO, QUINDI E' L'OSSERVABILE ESIGENTE
void estremi(Parametri par, std::vector<double> &psi,
             double &psiMin, double &psiMax)
{
    int M = par.N*par.N;
    psiMin = psi[0];
    psiMax = psi[0];
    for(int k=0;k<M;k++)
    {
        psiMin = std::min(psiMin,psi[k]);
        psiMax = std::max(psiMax,psi[k]);
    }
    psiMin = psiMin/(par.U*par.L);
    psiMax = psiMax/(par.U*par.L);
}


//=========================  SALVATAGGIO  =========================

//SALVATAGGIO DEL CAMPO PER IL PLOT IN PYTHON
//UNA RIGA PER CELLA, COORDINATE DEL CENTRO CELLA
void salvaCampo(Parametri par, std::vector<double> &u, std::vector<double> &v,
                std::vector<double> &p, std::vector<double> &psi,
                Risultato ris, double psiMin, double psiMax)
{
    int N = par.N;
    double dx = passo(par);

    FILE *fcampo = fopen("campo.csv","w");
    fprintf(fcampo,"x,y,u,v,p,psi\n");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            int k = j + N*i;
            fprintf(fcampo,"%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    (j+0.5)*dx,(i+0.5)*dx,u[k],v[k],p[k],psi[k]);
        }
    }
    fclose(fcampo);

    //PARAMETRI DELLA SIMULAZIONE, COSI' PYTHON NON DEVE INDOVINARLI
    FILE *fpar = fopen("parametri.csv","w");
    fprintf(fpar,"N,L,U,Re,rho,mu,gamma,au,ap,nswUV,nswP,iter,residuo,tempo,psiMin,psiMax\n");
    fprintf(fpar,"%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%.9g,%.9g,%.9g,%.9g\n",
            par.N,par.L,par.U,par.Re,par.rho,viscosita(par),par.gamma,
            par.au,par.ap,par.nswUV,par.nswP,
            ris.iter,ris.err,ris.tempo,psiMin,psiMax);
    fclose(fpar);

    printf("CAMPO SALVATO SU campo.csv E parametri.csv\n");
}


//====================  LETTURA DA RIGA DI COMANDO  ====================

//LEGGE ARGOMENTI DELLA FORMA chiave=valore.
//sscanf PROVA A LEGGERE DALLA STRINGA IL FORMATO INDICATO: SE CI RIESCE
//SCRIVE NELLA VARIABILE E RESTITUISCE 1, ALTRIMENTI RESTITUISCE 0.
//RESTITUISCE 1 SE E' STATO CHIESTO LO STUDIO SULLE SPAZZATE
int leggiArgomenti(int argc, char *argv[], Parametri &par)
{
    int studio = 0;

    for(int a=1;a<argc;a++)
    {
        if(sscanf(argv[a],"N=%d",&par.N)==1)           { continue; }
        if(sscanf(argv[a],"Re=%lf",&par.Re)==1)        { continue; }
        if(sscanf(argv[a],"U=%lf",&par.U)==1)          { continue; }
        if(sscanf(argv[a],"rho=%lf",&par.rho)==1)      { continue; }
        if(sscanf(argv[a],"au=%lf",&par.au)==1)        { continue; }
        if(sscanf(argv[a],"ap=%lf",&par.ap)==1)        { continue; }
        if(sscanf(argv[a],"gamma=%lf",&par.gamma)==1)  { continue; }
        if(sscanf(argv[a],"nswUV=%d",&par.nswUV)==1)   { continue; }
        if(sscanf(argv[a],"nswP=%d",&par.nswP)==1)     { continue; }
        if(sscanf(argv[a],"tol=%lf",&par.tol)==1)      { continue; }
        if(sscanf(argv[a],"maxIter=%d",&par.maxIter)==1) { continue; }
        if(sscanf(argv[a],"verbose=%d",&par.verbose)==1) { continue; }

        //sscanf CON UN FORMATO SENZA %  RESTITUISCE 0 SE LA STRINGA COMBACIA,
        //QUINDI PER LE OPZIONI SENZA VALORE SI CONFRONTA IL PRIMO CARATTERE
        if(argv[a][0]=='-' && argv[a][1]=='-' && argv[a][2]=='s')
        {
            studio = 1;
            continue;
        }

        printf("ARGOMENTO NON RICONOSCIUTO: %s\n",argv[a]);
        printf("USO: cavita [N=128] [Re=1000] [U=5] [rho=1000] [au=0.8] [ap=0.2]\n");
        printf("            [gamma=1] [nswUV=3] [nswP=20] [tol=1e-4] [maxIter=5000]\n");
        printf("            [verbose=1] [--spazzate]\n");
        printf("gamma: 1 = CDS, 0 = UDS, in mezzo una miscela\n");
    }

    return studio;
}


//==================  STUDIO SUL NUMERO DI SPAZZATE  ==================

//QUANTE SPAZZATE SIP CONVIENE SPENDERE SU OGNI SISTEMA A ITERAZIONE ESTERNA.
//NON HA SENSO RISOLVERE BENE I SISTEMI INTERNI, PERCHE' COEFFICENTI E TERMINE
//NOTO CAMBIERANNO ANCORA MOLTE VOLTE PRIMA CHE IL PROBLEMA NON LINEARE SIA
//RISOLTO (FERZIGER 12.2.2.1); MA FERMARSI TROPPO PRESTO FA CRESCERE IL NUMERO
//DI ITERAZIONI ESTERNE. L'OTTIMO DIPENDE DAL PROBLEMA E VA MISURATO.
//SI GUARDANO SIA LE ITERAZIONI ESTERNE SIA IL TEMPO: DIVERGONO, ED E'
//ESATTAMENTE QUELLO IL COMPROMESSO.
//SI PROVA SU DUE GRIGLIE, par.N/2 E par.N, PERCHE' L'OTTIMO SI SPOSTA CON
//IL RAFFINAMENTO: IL SISTEMA DI POISSON DIVENTA RELATIVAMENTE PIU' DURO
//MENTRE QUELLO DELLA QUANTITA' DI MOTO NO
void studioSpazzate(Parametri base)
{
    int listaP[6] = {2,5,10,20,30,50};
    int listaUV[5] = {1,2,3,5,10};
    int griglie[2];
    griglie[0] = base.N/2;
    griglie[1] = base.N;

    FILE *f = fopen("spazzate.csv","w");
    fprintf(f,"fase,N,nswUV,nswP,iter,tempo,convergiuto,res_u,res_v,res_massa\n");

    printf("STUDIO SULLE SPAZZATE, GRIGLIE %d e %d\n\n",griglie[0],griglie[1]);

    for(int g=0;g<2;g++)
    {
        Parametri par = base;
        par.N = griglie[g];
        par.verbose = 0;

        //FASE 1: SI VARIA nswP TENENDO FERMO nswUV
        printf("N = %d, fase 1: vario nswP con nswUV = %d\n",par.N,base.nswUV);
        printf("  nswP    iter    tempo(s)\n");

        int miglioreP = base.nswP;
        double migliorTempo = -1.0;

        for(int a=0;a<6;a++)
        {
            par.nswUV = base.nswUV;
            par.nswP = listaP[a];

            int M = par.N*par.N;
            std::vector<double> u(M,0.0), v(M,0.0), p(M,0.0);
            Risultato ris = risolviCavita(par,u,v,p);

            printf("  %4d   %5d    %7.2f%s\n",par.nswP,ris.iter,ris.tempo,
                   ris.convergiuto==1 ? "" : "   (NON CONVERGE)");
            fprintf(f,"nswP,%d,%d,%d,%d,%.4f,%d,%.6e,%.6e,%.6e\n",
                    par.N,par.nswUV,par.nswP,ris.iter,ris.tempo,ris.convergiuto,
                    ris.erru,ris.errv,ris.errm);

            if(ris.convergiuto==1 && (migliorTempo<0.0 || ris.tempo<migliorTempo))
            {
                migliorTempo = ris.tempo;
                miglioreP = par.nswP;
            }
        }

        //FASE 2: SI VARIA nswUV AL MIGLIOR nswP TROVATO
        printf("\nN = %d, fase 2: vario nswUV con nswP = %d\n",par.N,miglioreP);
        printf("  nswUV   iter    tempo(s)\n");

        for(int a=0;a<5;a++)
        {
            par.nswUV = listaUV[a];
            par.nswP = miglioreP;

            int M = par.N*par.N;
            std::vector<double> u(M,0.0), v(M,0.0), p(M,0.0);
            Risultato ris = risolviCavita(par,u,v,p);

            printf("  %4d   %5d    %7.2f%s\n",par.nswUV,ris.iter,ris.tempo,
                   ris.convergiuto==1 ? "" : "   (NON CONVERGE)");
            fprintf(f,"nswUV,%d,%d,%d,%d,%.4f,%d,%.6e,%.6e,%.6e\n",
                    par.N,par.nswUV,par.nswP,ris.iter,ris.tempo,ris.convergiuto,
                    ris.erru,ris.errv,ris.errm);
        }
        printf("\n");
    }

    fclose(f);
    printf("RISULTATI SALVATI SU spazzate.csv\n");
}


//============================  MAIN  ============================

int main(int argc, char *argv[])
{
    Parametri par = parametriDefault();
    int studio = leggiArgomenti(argc,argv,par);

    if(studio==1)
    {
        studioSpazzate(par);
        return 0;
    }

    int M = par.N*par.N;
    std::vector<double> u(M,0.0), v(M,0.0), p(M,0.0);

    printf("N = %d, Re = %g, gamma = %g (1 = CDS, 0 = UDS), nswUV = %d, nswP = %d\n\n",
           par.N,par.Re,par.gamma,par.nswUV,par.nswP);

    Risultato ris = risolviCavita(par,u,v,p);

    printf("\nTERMINATO DOPO %d ITERAZIONI, RESIDUO %10.3e, TEMPO %.2f s\n",
           ris.iter,ris.err,ris.tempo);
    if(ris.convergiuto==0)
    {
        printf("ATTENZIONE: USCITO PER maxIter, NON PER CONVERGENZA\n");
    }

    //INTENSITA' DEI VORTICI
    std::vector<double> om(M,0.0), psi(M,0.0), psiU(M,0.0), psiV(M,0.0);
    vorticita(par,u,v,om);
    double resPsi = funzioneCorrente(par,om,psi);

    double psiMin, psiMax;
    estremi(par,psi,psiMin,psiMax);

    //CONTROLLO: LE DUE INTEGRAZIONI DIRETTE DEVONO DARE LO STESSO psiMin
    double minU, maxU, minV, maxV;
    correnteIntegrata(par,u,v,0,psiU);
    correnteIntegrata(par,u,v,1,psiV);
    estremi(par,psiU,minU,maxU);
    estremi(par,psiV,minV,maxV);

    printf("\nINTENSITA' DEI VORTICI (psi/(U*L))\n");
    printf("  vortice primario     psi_min = %9.5f   (Ferziger 8.4.1: -0.11893)\n",psiMin);
    printf("  vortici secondari    psi_max = %9.5f   (Ferziger 8.4.1:  0.00173)\n",psiMax);
    printf("  controllo psi_min integrando u lungo y:  %9.5f\n",minU);
    printf("  controllo psi_min integrando v lungo x:  %9.5f\n",minV);
    printf("  residuo della Poisson per psi: %.1e\n",resPsi);
    if(resPsi>1e-8)
    {
        printf("  ATTENZIONE: POISSON NON RISOLTA, psi NON E' AFFIDABILE\n");
    }

    salvaCampo(par,u,v,p,psi,ris,psiMin,psiMax);

    return 0;
}
