#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include "CoreStructure.h"
#include "Tools.h"

namespace {
    void FCurrentinsert(map<vector<int>, int>* FCurrent, vector<int> currentxy, int* currentpos, string insertmode, vector<double>* symaxis) {
        if ( insertmode == "None" ) {
            ( *FCurrent ).insert(pair<vector<int>, int>(currentxy, *currentpos));
            ( *currentpos ) += 1;
        }
        else if ( insertmode == "4fold" ) {
            vector<int> sym1{ int(round(2 * ( *symaxis )[ 0 ] - currentxy[ 0 ])), currentxy[ 1 ] };
            vector<int> sym2{ int(round(2 * ( *symaxis )[ 0 ] - currentxy[ 0 ])), int(round(2 * ( *symaxis )[ 1 ] - currentxy[ 1 ])) };
            vector<int> sym3{ currentxy[ 0 ], int(round(2 * ( *symaxis )[ 1 ] - currentxy[ 1 ])) };
            if ( ( *FCurrent ).count(sym1) ) {
                ( *FCurrent ).insert(pair<vector<int>, int>(currentxy, ( *FCurrent )[ sym1 ]));
            }
            else if ( ( *FCurrent ).count(sym2) ) {
                ( *FCurrent ).insert(pair<vector<int>, int>(currentxy, ( *FCurrent )[ sym2 ]));
            }
            else if ( ( *FCurrent ).count(sym3) ) {
                ( *FCurrent ).insert(pair<vector<int>, int>(currentxy, ( *FCurrent )[ sym3 ]));
            }
            else {
                ( *FCurrent ).insert(pair<vector<int>, int>(currentxy, *currentpos));
                ( *currentpos ) += 1;
            }
        }
        else {
            cout << "FCurrentinsert: This sym mode not supported yet" << endl;
            throw 1;
        }
    }

    double calweight(int xo, int yo, int x, int y, double r) {
        return r - sqrt(double(x - xo) * double(x - xo) + double(y - yo) * double(y - yo));
    }

    bool circlerange(int xo, int yo, int x, int y, double r) {
        if ( double(x - xo) * double(x - xo) + double(y - yo) * double(y - yo) <= r * r - 0.01 ) {
            return true;
        }
        else {
            return false;
        }
    }

}  // namespace

CoreStructure::CoreStructure(double d_, VectorXi* geometry_, int Nx_, int Ny_, int Nz_, int N_, VectorXd* Inputdiel, bool Filter_, FilterOption* Filterstats_, string symmetry, vector<double> symaxis, bool Periodic_, int Lx_, int Ly_) {
    d = d_;
    Filter = Filter_;
    geometry = *geometry_;
    Nx = Nx_;
    Ny = Ny_;
    Nz = Nz_;
    N = N_;
    Periodic = Periodic_;
    Lx = Lx_;
    Ly = Ly_;

    int dividesym;
    if ( symmetry == "None" ) {
        dividesym = 1;
    }
    else if ( symmetry == "4fold" ) {
        dividesym = 4;
    }
    else {
        cout << "SpacePara: not None nor 4 fold. Not supported" << endl;
        throw 1;
    }

    int NFpara = int(round(( int(geometry.size( )) / 3 / Nz / dividesym )));    // number of free parameters. for extruded, symmetric, take one quadrant of one xy-plane of geo. 121 in standard case
    cout << "NFpara" << NFpara << endl;

    Para = VectorXd::Zero(NFpara);
    geometryPara = VectorXi::Zero(N);

    // stores an index that can be used to find the free parameter index associated with that index.
    // uses symmetry and reflections in FCurrentInsert to keep the range [0,120]. for example, 
    // geometryPara(Nx-1) = 1, geometry(Nx) = 0, geometry(Nx*Ny+1)=0 because of symmetry/extrusions.
    geometryPara = VectorXi::Zero(N);
    map<vector<int>, int> FCurrent;         // only used to design geometryPara, so can be removed if geometryPara designed differently
    int currentpos = 0;
    for ( int i = 0; i <= N - 1; i++ ) {
        int x = ( geometry ) ( 3 * i );
        int y = ( geometry ) ( 3 * i + 1 );

        vector<int> currentxy{ x,y };

        // for extruded geometries, same xy position can be used to refer to pixels that have different z-coord.
        if ( !FCurrent.count(currentxy) ) {

            FCurrentinsert(&FCurrent, currentxy, &currentpos, symmetry, &symaxis);
        }
        geometryPara(i) = FCurrent[ currentxy ];
        //cout << "geometryPara(" << i << ") is: " << geometryPara(i) << endl;
    }

    // This is used to map a free parameter position (0 to NFPara, or 121) to a vector of
    // positions that this free parameter maps to, considering relfections and extrusions.
    // for example, the first vector would be (0, Nx, Nx*Ny-Nx, Nx*Ny, ...). if Nz is 10,
    // and you have symmetry, each position will map to 10*4, or 40 other pixels.
    Paratogeometry = vector<vector<int>>(NFpara);
    for ( int i = 0; i <= N - 1; i++ ) {
        //cout << "dipole position : " << i << endl;

        ( Paratogeometry[ geometryPara(i) ] ).push_back(i);

        /*vector<int> currentvector = Paratogeometry[geometryPara(i)];
        for (int j = 0; j < currentvector.size(); j++) {
            cout << "Paratogeometry at index " << j << " is: " << currentvector[j] << endl;
        } */

    }

    // used to map a pixel vector to its inputdiel value (0-1)
    map<vector<int>, double> Inputmap;
    if ( geometry.size( ) != ( *Inputdiel ).size( ) ) {
        cout << "ERROR: SpacePara::SpacePara: Filter==(*InputGeo).size() != (*Inputdiel).size()" << endl;
        throw 1;
    }
    int Inputsize = int(round(int(( geometry ).size( )) / 3));
    for ( int i = 0; i < Inputsize; i++ ) {
        Inputmap.insert(pair<vector<int>, double>(vector<int>{( geometry ) ( 3 * i ), ( geometry ) ( 3 * i + 1 ), ( geometry ) ( 3 * i + 2 )}, ( *Inputdiel )( 3 * i )));
    }

    // used to fill Para, which are the parameter values of the free indices (NFpara)
    // or one quadrant in a symmetric, extruded structure. using Paratogeometry here
    // to fetch each position, but if this is the only use of Paratogeometry, seems useless.
    for ( int i = 0; i < NFpara; i++ ) {

        int pos = Paratogeometry[ i ][ 0 ];
        //cout << "Pos at position " << i << " is: " << pos << endl;
        vector<int> node{ ( geometry ) ( 3 * pos ), ( geometry ) ( 3 * pos + 1 ), ( geometry ) ( 3 * pos + 2 ) };
        Para(i) = Inputmap[ node ];
    }
    cout << "para values:" << endl;

    for ( int i = 0; i < Para.size( ); i++ ) {
        cout << Para(i) << " ";
    }

    if ( Filter == true ) {
        Para_origin = Para;
        Para_filtered = Para;
        Filterstats = Filterstats_;
        if ( Filterstats == NULL ) {
            cout << "ERROR: SpacePara::SpacePara: Filter==true then Filterstats must be passed in." << endl;
            throw 1;
        }

        assignFreeWeightsForFilter( );

    }

    // original CoreStructure stuff below

    //---------------------------------------------------initial diel------------------------------------
    diel_old = VectorXd::Zero(3 * N);
    diel_old_max = diel_old;
    for ( int i = 0; i <= N - 1; i++ ) {
        double dieltmp = Para(geometryPara(i));
        diel_old(3 * i) = dieltmp;
        diel_old(3 * i + 1) = dieltmp;
        diel_old(3 * i + 2) = dieltmp;
    }
}

void CoreStructure::UpdateStr(VectorXd step, int current_it, int Max_it) {
    cout << "step in UpdateStr: " << step.mean() << endl;

    int Parasize = Para.size();
    if (Parasize != step.size()) {
        cout << "ERROR: In CoreStructure::UpdateStr(VectorXd step), step.size!=FreePara.size";
        throw 1;
    }

    if (Filter == true) {
        //When there is filter
        (*Filterstats).update_beta(current_it, Max_it);                  //Update beta value according to current iteration
        
        for (int i = 0; i <= Parasize - 1; i++) {
            Para_origin(i) += step(i);
            if (Para_origin(i) >= 1) {
                Para_origin(i) = 1;
            }
            if (Para_origin(i) <= 0) {
                Para_origin(i) = 0;
            }
        }
        
        cout << "Beta at iteration " << current_it << " is " << Filterstats->get_beta() << endl;

        for (int i = 0; i <= Parasize - 1; i++) {
            int weightnum = (FreeWeight[i]).size();
            double numerator = 0.0;
            double denominator = 0.0;
            for (int j = 0; j <= weightnum - 1; j++) {
                numerator += (FreeWeight[i][j].weight) * Para_origin(FreeWeight[i][j].position);
                denominator += (FreeWeight[i][j].weight);

            }
            Para_filtered(i) = numerator / denominator;

            double Para_physical = Filterstats->SmoothDensity(Para_filtered(i));     
            Para(i) = Para_physical;

        }
    }
    else {//When there is no filter
        for (int i = 0; i <= Parasize - 1; i++) {
            Para(i) += step(i);
            if (Para(i) >= 1) {
                Para(i) = 1;
            }
            if (Para(i) <= 0) {
                Para(i) = 0;
            }
        }
    }


    for (int i = 0; i <= N - 1; i++) {
        int position = geometryPara(i);
        double value = Para(position);
        diel_old(3 * i) = value;
        diel_old(3 * i + 1) = value;
        diel_old(3 * i + 2) = value;
    }
}


void CoreStructure::UpdateStrSingle(int idx, double value) {

    diel_old(3 * idx) = value;
    diel_old(3 * idx + 1) = value;
    diel_old(3 * idx + 2) = value;

}

void CoreStructure::output_to_file(string save_position, int iteration, string mode) {

    if (mode == "normal") {
        string name;
        name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
        ofstream fout(name);
        fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
        fout << &geometry << endl;
        fout << diel_old << endl;
        fout << d << endl;
        fout.close();
    }
    else {
        string name;
        name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
        ofstream fout(name);
        fout << diel_old << endl;
        fout.close();
    }
}

void CoreStructure::assignFreeWeightsForFilter( ) {
    //Only works when freepara is 2D binding (Only do the filter in 2D)
    int NFpara = Para.size( );
    FreeWeight = vector<vector<WeightPara>>(NFpara);

    double rfilter = ( *Filterstats ).get_rfilter( );
    for ( int i = 0; i <= NFpara - 1; i++ ) {
        int poso = Paratogeometry[ i ][ 0 ];                 //As 2D extrusion is assumed, different z does not matter
        int xo = ( geometry ) ( 3 * poso );
        int yo = ( geometry ) ( 3 * poso + 1 );
        int zo = ( geometry ) ( 3 * poso + 2 );



        for ( int j = 0; j <= NFpara - 1; j++ ) {
            //bool inornot = false;
            //cout << j << endl;
            for ( int k = 0; k < Paratogeometry[ j ].size( ); k++ ) {
                int posr = Paratogeometry[ j ][ k ];

                int xr = ( geometry ) ( 3 * posr );
                int yr = ( geometry ) ( 3 * posr + 1 );
                int zr = ( geometry ) ( 3 * posr + 2 );
                if ( Periodic == false ) {
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr, yr, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                }
                else {
                    //Own cell
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr, yr, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //0, -1
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr, yr - Ly, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr - Ly, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //-1, -1
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr - Lx, yr - Ly, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr - Lx, yr - Ly, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //-1, 0
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr - Lx, yr, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr - Lx, yr, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //-1, 1
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr - Lx, yr + Ly, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr - Lx, yr + Ly, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //0, 1
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr, yr + Ly, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr, yr + Ly, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //1, 1
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr + Lx, yr + Ly, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr + Lx, yr + Ly, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //1, 0
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr + Lx, yr, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr + Lx, yr, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                    //1, -1
                    if ( ( zo == zr ) && ( circlerange(xo, yo, xr + Lx, yr - Ly, rfilter) ) ) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = j;
                        double weight = calweight(xo, yo, xr + Lx, yr - Ly, rfilter);
                        FreeWeight[ i ].push_back(WeightPara{ weight,posweight });
                        //break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j. But when there is symmetry, this is needed because same z can have differnt x, y.
                    }
                }


            }

        }
    }
}

VectorXi* CoreStructure::get_geometryPara( ) {
    return &geometryPara;
}

VectorXd* CoreStructure::get_Para( ) {
    return &Para;
}

VectorXd* CoreStructure::get_Para_origin( ) {
    if ( !Filter ) {
        cout << "ERROR: SpacePara::get_Para_origin()--Filter can not be false" << endl;
        throw 1;
    }
    return &Para_origin;
}

VectorXd* CoreStructure::get_Para_filtered( ) {
    if ( !Filter ) {
        cout << "ERROR: SpacePara::get_Para_filtered()--Filter can not be false" << endl;
        throw 1;
    }
    return &Para_filtered;
}

bool CoreStructure::get_Filter( ) {
    return Filter;
}

FilterOption* CoreStructure::get_Filterstats( ) {
    if ( !Filter ) {
        cout << "ERROR: SpacePara::get_Filterstats()--Filter can not be false" << endl;
        throw 1;
    }
    return Filterstats;
}

vector<vector<WeightPara>>* CoreStructure::get_FreeWeight( ) {
    if ( !Filter ) {
        cout << "ERROR: SpacePara::get_FreeWeight()---Filter can not be false" << endl;
        throw 1;
    }
    return &FreeWeight;
}

vector<vector<int>>* CoreStructure::get_Paratogeometry( ) {
    return &Paratogeometry;
}

// from StructureAndSpace.cpp

VectorXd* CoreStructure::get_diel_old() {
    return &diel_old;
}
VectorXd* CoreStructure::get_diel_old_max() {
    return &diel_old_max;
}

