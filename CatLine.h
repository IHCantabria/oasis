/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CLASE CatLine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
class CatLine {
public:
	int nLine, nNodos;
	double L, dL, EA, beta, rho0, d, A, Cdt, Cdn, Cmn, CB, GK, GC, Gmu, Gvc, Dz;
	double tenAnch[3], tenFair[3], posAnch[3], posFair[3];
	double * pos, * vel, * acc;
	double xF, zF, HF, VF;

	void set_nLine(int n){nLine=n;}
	void leer_datosMoorings(void);
	void initLine(void);
	void qs_GetTen(void);
};

