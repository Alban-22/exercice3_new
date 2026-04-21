#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>
#include <valarray>        // input output manipulators
#include "common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>
#include <utility> // pour std::pairf

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private:
    // Existing private members of Engine...
  const double pi=3.1415926535897932384626433832795028841971e0;
  const double G=6.67430e-11; // Constante de gravitation universelle
  // definition des variables
  unsigned int N;
  bool adaptive;
  unsigned int  n; 
  double t;  // Temps courant pas de temps
  double tf;          // Temps final
  double dt;
  double dt_current;      // Intervalle de temps
  int N_excit;  // Nombre de périodes d'excitation
  int nsteps_per; // Nombre de pas de temps par période d'excitation
  double epsilon;
  double s;
  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie
  double L, r, kappa, d;         // accélération gravitationnelle, masse, longueur, fréquence angulaire, rayon, coefficient de frottement
  valarray<double> y;
  valarray<double> S;
  valarray<double> R;
  valarray<double> m;    //k1 de RK4
  double rho_0;
  double R_T;
  double lambda;
  valarray<double> Cx;
  valarray<double> alpha;
  valarray<bool> bool_alpha;
  valarray<double> vnorm;
  

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  // TODO definir l'énergie mecanique
  

  double ix(int i) const {
    return y[4*(i)];
  }
  double iy(int i) const {
    return y[4*(i) +1];
  }
  double ivx(int i) const {
    return y[4*(i) +2];
  }
  double ivy(int i) const {
    return y[4*(i) +3];
  }

  // TODO definir la puissance des forces non conservatives
double compute_norm(double x, double y){ return sqrt(x*x+y*y);};


std::pair<double,double> compute_F_ij_grav(size_t i, size_t j)
{
    double dx = ix(i) - ix(j);
    double dy = iy(i) - iy(j);
    double dist = compute_norm(dx, dy);

    if (dist < 1e-12) {
        return {0.0, 0.0};
    }

    double factor = -G * m[i] * m[j] / (dist*dist*dist);

    return {factor * dx, factor * dy};
}

double compute_rho(double r){ 
  if (r<R_T) { r=R_T;}
  return rho_0 * exp(-(r-R_T)/lambda); 
};


std::pair<double,double> compute_F_ij_frott(size_t i, size_t j)
{
    if (Cx[j] == 0 || S[i] == 0) {
        return {0.0, 0.0};
    }

    double dvx = ivx(i) - ivx(j);
    double dvy = ivy(i) - ivy(j);
    double vrel = compute_norm(dvx, dvy);

    if (vrel < 1e-12) {
        return {0.0, 0.0};
    }

    double rho = compute_rho(compute_norm(ix(i)-ix(j), iy(i)-iy(j)));

    double factor = -0.5 * rho * S[i] * Cx[j] * vrel;

    return {factor * dvx, factor * dvy};
}


  // TODO écrire la fonction pour l'acceleration

  
std::pair<double,double> compute_acc_sur_i(size_t i){
    double Fx_tot = 0.0;
    double Fy_tot = 0.0;

    for (size_t j = 0; j < N + n; ++j) {
        if (j == i) continue;

        auto [Fx_grav, Fy_grav]   = compute_F_ij_grav(i, j);
        auto [Fx_frott, Fy_frott] = compute_F_ij_frott(i, j);

        Fx_tot += Fx_grav + Fx_frott;
        Fy_tot += Fy_grav + Fy_frott;
    }

    return {Fx_tot / m[i], Fy_tot / m[i]};
}

double compute_acc_sur_i_norm(size_t i) {
    auto [ax, ay] = compute_acc_sur_i(i);
    return sqrt(ax*ax + ay*ay);
}

double compute_energie_cinetique() const
{
    double Ec = 0.0;
    for (size_t i = 0; i < N + n; ++i) {
        double vx = ivx(i);
        double vy = ivy(i);
        Ec += 0.5 * m[i] * (vx*vx + vy*vy);
    }
    return Ec;
}


double compute_energie_potentielle_grav() const
{
    double Ep = 0.0;
    for (size_t i = 0; i < N + n; ++i) {
        for (size_t j = i + 1; j < N + n; ++j) {
            double dx = ix(i) - ix(j);
            double dy = iy(i) - iy(j);
            double rij = sqrt(dx*dx + dy*dy);

            if (rij > 1e-12) {
                Ep += -G * m[i] * m[j] / rij;
            }
        }
    }
    return Ep;
}


double compute_energie_mecanique() const
{
    return compute_energie_cinetique() + compute_energie_potentielle_grav();
}

std::pair<double,double> compute_quantite_mouvement() const
{
    double Px = 0.0;
    double Py = 0.0;

    for (size_t i = 0; i < N + n; ++i) {
        Px += m[i] * ivx(i);
        Py += m[i] * ivy(i);
    }

    return {Px, Py};
}

double compute_quantite_mouvement_norm() const{

    auto [Px, Py] = compute_quantite_mouvement();
    return sqrt(Px*Px + Py*Py);
}


bool checkCollisions() const
{
    for (size_t i = 0; i < N + n; ++i) {
        for (size_t j = i + 1; j < N + n; ++j) {

            double dx = ix(i) - ix(j);
            double dy = iy(i) - iy(j);
            double dist = sqrt(dx*dx + dy*dy);

            if (dist <= R[i] + R[j]) {
                return true;
            }
        }
    }
    return false;
}

 void compute_f(valarray<double>& f){
        for (unsigned int i(0); i<N+n; ++i){
             f[4*i] = ivx(i); // vx 
             f[4*i+1] = ivy(i); // vy
             auto [ax, ay] = compute_acc_sur_i(i);
             f[4*i+2] = ax; // ax
             f[4*i+3] = ay; // ay
             
        }
  } 

valarray<double> step_fixe(double dt_local, valarray<double> y_in)
{
    // Sauvegarde de l'état courant de la classe
    valarray<double> y_save = y;

    // Vecteurs RK4
    valarray<double> k1(0.0, y_in.size());
    valarray<double> k2(0.0, y_in.size());
    valarray<double> k3(0.0, y_in.size());
    valarray<double> k4(0.0, y_in.size());


    
    // k1 = f(y_in)
    y = y_in;
    compute_f(k1);

    // k2 = f(y_in + dt/2 * k1)
    y = y_in + 0.5 * dt_local * k1;
    compute_f(k2);

    // k3 = f(y_in + dt/2 * k2)
    y = y_in + 0.5 * dt_local * k2;
    compute_f(k3);

    // k4 = f(y_in + dt * k3)
    y = y_in + dt_local * k3;
    compute_f(k4);

    // Résultat RK4
    valarray<double> y_out = y_in + (dt_local / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    // Restauration de l'état original de la classe
    y = y_save;

    return y_out;
}
  
void step() {
    y = step_fixe(dt_current, y);
    t += dt_current;
}
void step_adaptatif()
{
    double dt_try = dt_current;

    while (true) {
        valarray<double> y_half  = step_fixe(0.5 * dt_try, y);
        valarray<double> y_2half = step_fixe(0.5 * dt_try, y_half);
        valarray<double> y_full  = step_fixe(dt_try, y);

        valarray<double> delta_y = (y_full - y_2half);
        double err2= (delta_y * delta_y).sum();
        double err = sqrt(err2) / 15.0;

        if (err <= epsilon || err < 1e-30) {
            y = y_2half;     // meilleure approximation
            t += dt_try;

            if (err < 1e-30) {
                dt_current = 2.0 * dt_try;
            } else {
                dt_current = s * dt_try * pow(epsilon / err, 0.2);
            }
            break;
        }

        dt_try = s * dt_try * pow(epsilon / err, 0.2);
    }
}
  

void printOut(bool write)
{
  // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
  if((!write && last >= sampling) || (write && last != 1))
  {
    *outputFile << setprecision(15) << fixed
                << t << " ";

    for (size_t i = 0; i < N + n; ++i)
    {
      *outputFile << setprecision(15) << fixed
                  << ix(i) << " " << iy(i) << " "
                  << ivx(i) << " " << ivy(i) << " "
                  << compute_acc_sur_i_norm(i) << " ";
    }

    *outputFile << setprecision(15) << fixed
                << compute_energie_cinetique() << " "
                << compute_energie_potentielle_grav() << " "
                << compute_energie_mecanique() << " "
                << compute_quantite_mouvement_norm() << " "
                << dt_current
                << endl;

    last = 1;
  }
  else
  {
    last++;
  }
}


public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Read parameters from the config file
      N = configFile.get<unsigned int>("N");
      n = configFile.get<unsigned int>("n");
      tf = configFile.get<double>("tf");
      dt = configFile.get<double>("dt");
      dt_current = dt;
      N_excit = configFile.get<int>("N_excit");
      nsteps_per = configFile.get<int>("nsteps_per");
      sampling = configFile.get<unsigned int>("sampling");
      adaptive = configFile.get<bool>("adaptive");
      L = configFile.get<double>("L");
      r = configFile.get<double>("r");
      kappa = configFile.get<double>("kappa");
      d = configFile.get<double>("d");
      rho_0 = configFile.get<double>("rho_0");
      R_T = configFile.get<double>("R_T");
      lambda = configFile.get<double>("lambda");
      epsilon = configFile.get<double>("epsilon");
      s       = configFile.get<double>("s");  
      y = valarray<double>(0.0, 4*(N+n));
      S = valarray<double>(0.0, N+n);
      m = valarray<double>(0.0, N+n);
      R = valarray<double>(0.0, N+n);
      Cx = valarray<double>(0.0, N+n);
      alpha= valarray<double>(0.0, N+n);
      bool_alpha= valarray<bool>(false,N+n);
      vnorm = valarray<double>(0.0, N+n);
      for (unsigned int i = 0; i < N+n; ++i) {
          m[i] = configFile.get<double>("m" + to_string(i));
          S[i] = configFile.get<double>("S" + to_string(i));
          R[i]  = configFile.get<double>("R"  + to_string(i));
          Cx[i] = configFile.get<double>("Cx" + to_string(i));
          alpha[i] = configFile.get<double>("alphadeg" + to_string(i),0.0) * pi / 180.0; // Convertir en radians
          bool_alpha[i] = configFile.get<bool>("bool_alpha" + to_string(i), false);
          vnorm[i] = configFile.get<double>("vnorm" + to_string(i), 0.0);
          y[4*i]   = configFile.get<double>("x" + to_string(i));
          y[4*i+1] = configFile.get<double>("y" + to_string(i));
          if(bool_alpha[i]) {
            y[4*i+2] = vnorm[i] * cos(alpha[i]);
            y[4*i+3] = vnorm[i] * sin(alpha[i]);
          }
          else {
            y[4*i+2] = configFile.get<double>("vx" + to_string(i));
            y[4*i+3] = configFile.get<double>("vy" + to_string(i));
          }
          
          
}

      // Initialize the output file
      string output_name = configFile.get<string>("output", "output.txt");
      outputFile = new ofstream(output_name);
    };


    // Destructeur virtuel
    virtual ~Engine()
    {
      outputFile->close();
      delete outputFile;
    };
      // Simulation complete

void run()
{
    t = 0.;
    last = 0;
    printOut(true);

    while (t < tf)
    {
        double dt_nominal = dt_current;
        if (t + dt_current > tf) {
            dt_current = tf - t;
        }

        if (adaptive) {
            step_adaptatif();
        } else {
            step();
        }

        printOut(false);

        if (checkCollisions()) {
            cout << "Collision detectee a t = " << t << " s  !!!!!!" << endl;
            break;
        }

        if (!adaptive) {
            dt_current = dt_nominal;
        }
    }

    printOut(true);
}
};
// programme
int main(int argc, char* argv[])
{
  // Existing main function implementation
  // ...
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

  Engine* engine;

  // Create an instance of Engine instead of EngineEuler
  engine = new Engine(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}
