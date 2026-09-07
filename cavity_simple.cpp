//CAVITA' QUADRATA LID-DRIVEN CON VELOCITA' COSTANTE U
#include<cstdio>
#include<cmath>
#include<vector>
#include<algorithm>

const double L = 1.0;
const double U = 5.0;
const double Re = 1000.0;
const double rho = 1000.0;
const double au = 0.8;
const double ap = 0.2;
const double tol = 1e-5;

//METODO FV CON ALGORITMO SIMPLE,GRIGLIA CARTESIANA UNIFORME,DISCRETIZZAZIONE CDS


//METODO SIP DI STONE PER MATRICI A 5 DIAGONALI

double SIP(int N,
           std::vector<double> &aP, std::vector<double> &aE, std::vector<double> &aW,
           std::vector<double> &aN, std::vector<double> &aS, std::vector<double> &q,
           std::vector<double> &f, int nsw)
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

            LS[k] = aS[k]/(1.0 + alfa*UEs);
            LW[k] = aW[k]/(1.0 + alfa*UNw);

            double p1 = alfa*LS[k]*UEs;
            double p2 = alfa*LW[k]*UNw;

            LP[k] = 1.0/(aP[k] + p1 + p2 - LS[k]*UNs - LW[k]*UEw);   //LP = 1/L^P

            UE[k] = (aE[k] - p1)*LP[k];
            UN[k] = (aN[k] - p2)*LP[k];
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

                double r = q[k] - aP[k]*f[k];
                if(j!=N-1) r = r - aE[k]*f[k+1];
                if(j!=0)   r = r - aW[k]*f[k-1];
                if(i!=N-1) r = r - aN[k]*f[k+N];
                if(i!=0)   r = r - aS[k]*f[k-N];

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


int main()
{
    int N = 128;
    int M = N*N;
    int iter = 0;
    double dx = L/N;
    double mu = rho*U*L/Re;
    double resu, resv, resm;        //RESIDUI GREZZI, resu E resv ARRIVANO DAL SIP
    double erru, errv, errm, err;   //GLI STESSI, NORMALIZZATI

    //CAMPI AL CENTRO CELLA, TUTTI N x N A ZERO
    std::vector<double> riga(N,0.0);
    std::vector<std::vector<double>> u(N,riga), v(N,riga), p(N,riga);
    std::vector<std::vector<double>> uc(N,riga), vc(N,riga);
    std::vector<std::vector<double>> us(N,riga), vs(N,riga), pc(N,riga);

    //COEFFICENTI: x, y, PRESSIONE
    std::vector<double> aE(M,0.0), aW(M,0.0), aN(M,0.0), aS(M,0.0), aP(M,0.0), qP(M,0.0);
    std::vector<double> aEy(M,0.0), aWy(M,0.0), aNy(M,0.0), aSy(M,0.0), aPy(M,0.0), qPy(M,0.0);
    std::vector<double> aEp(M,0.0), aWp(M,0.0), aNp(M,0.0), aSp(M,0.0), aPp(M,0.0), qPp(M,0.0);

    //PORTATE SULLE FACCE, LORO CORREZZIONI, INCOGNITE DEI TRE SISTEMI
    std::vector<double> me(M,0.0), mw(M,0.0), mn(M,0.0), ms(M,0.0), dmp(M,0.0);
    std::vector<double> ue(M,0.0), vn(M,0.0);
    std::vector<double> ustar(M,0.0), vstar(M,0.0), pstar(M,0.0);

    //PORTATE MEMORIZZATE SULLA FACCIA EST E NORD DI OGNI CELLA.
    //NON SI RICALCOLANO DA ZERO OGNI VOLTA: SONO UNO STATO CHE SOPRAVVIVE
    //TRA UN'ITERAZIONE E L'ALTRA, ALTRIMENTI LA CORREZZIONE DI p' VA PERSA.
    //gpx E gpy SONO I GRADIENTI DI p AL CENTRO CELLA, SERVONO A RHIE-CHOW
    std::vector<double> fe(M,0.0), fn(M,0.0), gpx(M,0.0), gpy(M,0.0);

    do
    {
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
              //COEFFICENTI PER LA COMPONENTE X DELL'EQUAZIONE
              int k = j + N*i;

              //LA PORTATA DI FACCIA NON SI INTERPOLA PIU' QUI: SI LEGGE DA fe/fn,
              //CHE SONO COSTRUITE CON RHIE-CHOW DOPO LA QUANTITA' DI MOTO.
              //LA FACCIA OVEST DELLA CELLA k E' LA FACCIA EST DELLA CELLA k-1
              if(j!=N-1) { aE[k] =  0.5*fe[k]   - mu; }
              if(j!=0)   { aW[k] = -0.5*fe[k-1] - mu; }
              if(i!=N-1) { aN[k] =  0.5*fn[k]   - mu; }
              if(i!=0)   { aS[k] = -0.5*fn[k-N] - mu; }

              //PRESSIONE SULLE FACCE EST E OVEST
              //SE LA CELLA E' A PARETE SI USA LA PRESSIONE DELLA CELLA STESSA
              double pW = p[i][j];
              double pE = p[i][j];
              if(j!=0)   { pW = p[i][j-1]; }
              if(j!=N-1) { pE = p[i][j+1]; }
              qP[k] = dx*0.5*(pW - pE);

              //NO-SLIP. LA PARETE STA A dx/2 DAL CENTRO CELLA, NON A dx,
              //QUINDI IL SUO COEFFICENTE DIFFUSIVO E' mu*dx/(dx/2) = 2*mu.
              //ATTRAVERSO LA PARETE NON PASSA MASSA, QUINDI NIENTE CONVEZIONE.
              //LA PARTE INCOGNITA DEL FLUSSO VA SU aP, QUELLA NOTA SU qP.
              double aWall = 0.0;
              if(j==N-1) { aWall = aWall + 2*mu; }
              if(j==0)   { aWall = aWall + 2*mu; }
              if(i==0)   { aWall = aWall + 2*mu; }
              if(i==N-1) { aWall = aWall + 2*mu; qP[k] = qP[k] + 2*mu*U; }  //COPERCHIO

              //PORTATA NETTA USCENTE DALLA CELLA, SERVE PER aP
              double dm = 0.0;
              if(j!=N-1) { dm = dm + fe[k];   }
              if(j!=0)   { dm = dm - fe[k-1]; }
              if(i!=N-1) { dm = dm + fn[k];   }
              if(i!=0)   { dm = dm - fn[k-N]; }

              //aP E' MENO LA SOMMA DEI VICINI PIU' LA PORTATA NETTA
              //I VICINI CHE NON ESISTONO VALGONO 0 E NON DANNO CONTRIBUTO
              aP[k] = -(aE[k] + aW[k] + aN[k] + aS[k]) + dm + aWall;

              //COEFFICENTI PER LA COMPONENTE Y DELL'EQUAZIONE
              //STESSE FACCE E STESSE PORTATE DELL'EQUAZIONE IN x
              if(j!=N-1) { aEy[k] =  0.5*fe[k]   - mu; }
              if(j!=0)   { aWy[k] = -0.5*fe[k-1] - mu; }
              if(i!=N-1) { aNy[k] =  0.5*fn[k]   - mu; }
              if(i!=0)   { aSy[k] = -0.5*fn[k-N] - mu; }

              //PRESSIONE SULLE FACCE NORD E SUD
              double pS = p[i][j];
              double pN = p[i][j];
              if(i!=0)   { pS = p[i-1][j]; }
              if(i!=N-1) { pN = p[i+1][j]; }
              qPy[k] = dx*0.5*(pS - pN);

              //dm E aWall SONO GLI STESSI DI PRIMA: DIPENDONO DAL CAMPO E DALLA
              //GEOMETRIA, NON DALL'EQUAZIONE. QUI v DI PARETE E' 0 OVUNQUE,
              //ANCHE SUL COPERCHIO, QUINDI qPy NON RICEVE NIENTE
              aPy[k] = -(aEy[k] + aWy[k] + aNy[k] + aSy[k]) + dm + aWall;

              //SOTTORILASSAMENTO
              aP[k] = aP[k]/au;
              aPy[k] = aPy[k]/au;
              qP[k] = qP[k] + (1-au)*aP[k]*u[i][j];
              qPy[k] = qPy[k] + (1-au)*aPy[k]*v[i][j];
            }
        }

        //RISOLUZIONE SISTEMA LINEARE CON METODO SIP
        //SI PARTE DAL CAMPO DELL'ITERAZIONE PRECEDENTE COME STIMA INIZIALE
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                ustar[k] = u[i][j];
                vstar[k] = v[i][j];
            }
        }
        //POCHE SPAZZATE: DENTRO AL SIMPLE NON SERVE RISOLVERE A CONVERGENZA
        resu = SIP(N,aP,aE,aW,aN,aS,qP,ustar,3);
        resv = SIP(N,aPy,aEy,aWy,aNy,aSy,qPy,vstar,3);

        //GRADIENTI DI p AL CENTRO CELLA (STESSE FACCE USATE PER qP e qPy)
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                double pW = p[i][j], pE = p[i][j], pS = p[i][j], pN = p[i][j];
                if(j!=0)   { pW = p[i][j-1]; }
                if(j!=N-1) { pE = p[i][j+1]; }
                if(i!=0)   { pS = p[i-1][j]; }
                if(i!=N-1) { pN = p[i+1][j]; }
                gpx[k] = (pE - pW)*0.5/dx;
                gpy[k] = (pN - pS)*0.5/dx;
            }
        }

        //PORTATE SULLE FACCE CON INTERPOLAZIONE DI QUANTITA' DI MOTO (RHIE-CHOW).
        //LA MEDIA (u_P+u_E)/2 SI PORTA DIETRO I GRADIENTI CENTRATI DELLE DUE CELLE,
        //NEI QUALI (p_E - p_P) SI SEMPLIFICA: LA PORTATA NON VEDREBBE IL SALTO DI
        //PRESSIONE ATTRAVERSO LA FACCIA CHE STA ATTRAVERSANDO. LA PARENTESI TOGLIE
        //IL GRADIENTE MEDIATO E RIMETTE QUELLO A DUE PUNTI. SU CAMPO REGOLARE I DUE
        //QUASI COINCIDONO (TERMINE O(dx^3)), SU UNA SCACCHIERA NO E LA SMORZANO.
        //COSI' LA PORTATA RISPONDE A (p_E - p_P) CON COEFFICENTE (V/aP)_e, CHE E'
        //ESATTAMENTE QUELLO SCRITTO IN aEp: LE DUE EQUAZIONI SONO COERENTI
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                if(j!=N-1)
                {
                double C = 0.5*dx*dx*(1.0/aP[k] + 1.0/aP[k+1]);
                fe[k] = rho*dx*( 0.5*(ustar[k] + ustar[k+1])
                      - C*((p[i][j+1] - p[i][j])/dx - 0.5*(gpx[k] + gpx[k+1])) );
                }
                if(i!=N-1)
                {
                double C = 0.5*dx*dx*(1.0/aPy[k] + 1.0/aPy[k+N]);
                fn[k] = rho*dx*( 0.5*(vstar[k] + vstar[k+N])
                      - C*((p[i+1][j] - p[i][j])/dx - 0.5*(gpy[k] + gpy[k+N])) );
                }
            }
        }

        //CALCOLO PORTATA MASSICA IN ECCEDENZA
        //LE FACCE DI PARETE RESTANO A 0 PERCHE' NON VENGONO MAI SCRITTE
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                if(j!=N-1) { me[k] = fe[k];   }
                if(j!=0)   { mw[k] = fe[k-1]; }
                if(i!=N-1) { mn[k] = fn[k];   }
                if(i!=0)   { ms[k] = fn[k-N]; }
                dmp[k] = mn[k] - ms[k] + me[k] - mw[k];
            }
        }

        //RISOLUZIONE EQUAZIONE DELLA CORREZZIONE DELLA PRESSIONE
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                //IL COEFFICENTE E' -rho*(V/aP) DELLA FACCIA, E LA FACCIA E' UNA SOLA:
                //SI USA LA MEDIA DELLE DUE CELLE, COSI' aEp[k] E aWp[k+1] COINCIDONO.
                //ATTRAVERSO UNA PARETE NON C'E' CORREZZIONE, QUINDI IL COEFFICENTE E' 0
                if(j!=N-1) { aEp[k] = -rho*0.5*dx*dx*(1.0/aP[k]  + 1.0/aP[k+1]);  }
                if(j!=0)   { aWp[k] = -rho*0.5*dx*dx*(1.0/aP[k]  + 1.0/aP[k-1]);  }
                if(i!=N-1) { aNp[k] = -rho*0.5*dx*dx*(1.0/aPy[k] + 1.0/aPy[k+N]); }
                if(i!=0)   { aSp[k] = -rho*0.5*dx*dx*(1.0/aPy[k] + 1.0/aPy[k-N]); }
                aPp[k] = -(aEp[k] + aWp[k] + aNp[k] + aSp[k]);
            }
        }
        //RISOLUZIONE SISTEMA LINEARE CON METODO SIP
        //IL TERMINE NOTO E' MENO LA PORTATA IN ECCEDENZA, E SI RIPARTE SEMPRE DA p'=0
        for(int k=0;k<M;k++)
        {
            qPp[k] = -dmp[k];
            pstar[k] = 0.0;
        }
        //LA PRESSIONE E' DEFINITA A MENO DI UNA COSTANTE, QUINDI LA MATRICE E'
        //SINGOLARE: SI FISSA p' = 0 IN UNA CELLA DI RIFERIMENTO (LA 0)
        aEp[0] = 0.0;
        aWp[0] = 0.0;
        aNp[0] = 0.0;
        aSp[0] = 0.0;
        aPp[0] = 1.0;
        qPp[0] = 0.0;

        //QUI SERVONO PIU' SPAZZATE: E' L'EQUAZIONE ELLITTICA
        SIP(N,aPp,aEp,aWp,aNp,aSp,qPp,pstar,20);

        //CALCOLO CORREZZIONI DELLE PORTATE DI FACCIA
        //m' = -rho*(V/aP)*(p'_E - p'_P), CHE E' GIA' IL COEFFICENTE aEp:
        //NON SERVE RICALCOLARLO. VANNO SOMMATE ALLE PORTATE MEMORIZZATE,
        //ED E' IL PASSO CHE RENDE dmp NULLO A CONVERGENZA
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                if(j!=N-1) { ue[k] = aEp[k]*(pstar[k+1] - pstar[k]); fe[k] = fe[k] + ue[k]; }
                if(i!=N-1) { vn[k] = aNp[k]*(pstar[k+N] - pstar[k]); fn[k] = fn[k] + vn[k]; }
            }
        }

        //CALCOLO CAMPO DELLE VELOCITA' CORRETTE
        //STESSA FORMULA DELLE FACCE MA COL GRADIENTE CENTRATO SULLA CELLA
        //SULLA PARETE p' DI FACCIA E' QUELLO DELLA CELLA, COME PER p
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                int k = j + N*i;
                double pcW = pstar[k];
                double pcE = pstar[k];
                double pcS = pstar[k];
                double pcN = pstar[k];
                if(j!=0)   { pcW = pstar[k-1]; }
                if(j!=N-1) { pcE = pstar[k+1]; }
                if(i!=0)   { pcS = pstar[k-N]; }
                if(i!=N-1) { pcN = pstar[k+N]; }

                uc[i][j] = dx*0.5*(pcW - pcE)/aP[k];
                vc[i][j] = dx*0.5*(pcS - pcN)/aPy[k];

                us[i][j] = ustar[k];
                vs[i][j] = vstar[k];
                pc[i][j] = pstar[k];
            }
        }

        //SOMMA TRA CAMPO DELLE VELOCITA' STAR E QUELLE CORRETTE
        //CALCOLO PRESSIONE CON SOTTORILASSAMENTO
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                u[i][j] = us[i][j] + uc[i][j];
                v[i][j] = vs[i][j] + vc[i][j];
                p[i][j] = p[i][j] + pc[i][j]*ap;
            }
        }

        //CRITERIO DI ARRESTO SUI RESIDUI
        
        //resu E resv LI RESTITUISCE GIA' IL SIP: SONO LA SOMMA DI |q - A*f|
        //CALCOLATA PRIMA DI RISOLVERE, CIOE' DI QUANTO L'EQUAZIONE DISCRETA
        //NON E' SODDISFATTA DAL CAMPO CORRENTE
        //resm E' IL RESIDUO DELLA CONTINUITA': LA SOMMA DEGLI SBILANCI DI MASSA
        resm = 0.0;
        for(int k=0;k<M;k++)
        {
            resm = resm + std::abs(dmp[k]);
        }

        //I RESIDUI SONO DIMENSIONALI, QUINDI VANNO NORMALIZZATI CON I VALORI
        //DI RIFERIMENTO DEL PROBLEMA, ALTRIMENTI tol DIPENDE DA rho, U, L
        erru = resu/(rho*U*U*L);
        errv = resv/(rho*U*U*L);
        errm = resm/(rho*U*L);

        //SI USA IL PIU' GRANDE DEI TRE
        err = std::max(erru,errv);
        err = std::max(err,errm);

        iter = iter + 1;
        printf("iter %5d   res_u %10.3e   res_v %10.3e   res_massa %10.3e\n",
               iter,erru,errv,errm);

        //TETTO ALLE ITERAZIONI PER NON RESTARE BLOCCATI SE LA SOLUZIONE DIVERGE
    }while(err>tol && iter<5000);

    printf("\nTERMINATO DOPO %d ITERAZIONI, RESIDUO %10.3e\n",iter,err);

    //SALVATAGGIO DEL CAMPO PER IL PLOT IN PYTHON
    //UNA RIGA PER CELLA, COORDINATE DEL CENTRO CELLA
    FILE *fcampo = fopen("campo.csv","w");
    fprintf(fcampo,"x,y,u,v,p\n");
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            fprintf(fcampo,"%.9g,%.9g,%.9g,%.9g,%.9g\n",
                    (j+0.5)*dx,(i+0.5)*dx,u[i][j],v[i][j],p[i][j]);
        }
    }
    fclose(fcampo);

    //PARAMETRI DELLA SIMULAZIONE, COSI' PYTHON NON DEVE INDOVINARLI
    FILE *fpar = fopen("parametri.csv","w");
    fprintf(fpar,"N,L,U,Re,rho,mu,iter,residuo\n");
    fprintf(fpar,"%d,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%.9g\n",N,L,U,Re,rho,mu,iter,err);
    fclose(fpar);

    printf("CAMPO SALVATO SU campo.csv E parametri.csv\n");
}
