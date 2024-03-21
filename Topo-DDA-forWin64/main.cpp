#define PI 3.14159265
#define _USE_MATH_DEFINES

#include <chrono>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <random>

#include "EvoDDAModel.h"
#include "filterReader.h"
#include "INIReader.h"
#include "symReader.h"
#include "ObjReader.h"
#include "Tools.h"

using namespace std::chrono;
namespace fs = std::filesystem;


ObjDDAModel* ObjFactoryOut(string ObjectName, vector<double> ObjectParameters, DDAModel* ObjDDAModel) {
    /*if (HavePenalty) {
        cout << "Using L1 Penalty with Penalty Factor " << PenaltyFactor << endl;
    }*/
    if (ObjectName == "PointE") {
        return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
    }
    if (ObjectName == "IntegratedE") {
        return new ObjIntegratedEDDAModel(ObjectParameters, ObjDDAModel);
    }
    cout << "NOT A LEGIT OBJECTIVE NAME!" << endl;
    return new ObjPointEDDAModel(ObjectParameters, ObjDDAModel);
}

void task() {
    INIReader reader1("task.ini");

    if (reader1.ParseError() != 0) {
        std::cout << "Can't load 'task.ini'\n";
        return;
    }

    map<string, string> findtask;
    findtask.insert(pair<string, string>("DDA verify", "DDA_verify_path"));
    findtask.insert(pair<string, string>("DDA input", "DDA_input_path"));
    findtask.insert(pair<string, string>("EvoOpt", "EvoOpt_path"));
    findtask.insert(pair<string, string>("EvoOpt 2D input", "EvoOpt_2D_input_path"));
    findtask.insert(pair<string, string>("EvoOpt 2D input periodic", "EvoOpt_2D_input_periodic_path"));
    findtask.insert(pair<string, string>("NN data generate", "NN_path"));

    string tasktype = reader1.Get("Model Name", "name", "UNKNOWN");
    string tasktypepath = findtask[tasktype];
    string path = reader1.Get("Path", tasktypepath, "UNKNOWN");

    cout << path << endl;
    INIReader reader2(path);

    // HEEYOO!!!
    if (tasktype == "EvoOpt 2D input periodic") {

        std::vector<double> epsilonArray = {0.5};
        std::vector<double> weightArray = {0};
        std::vector<int> projectionArray = { 5, 10, 20, 50 };
        std::vector<string> penaltytypeArray = {"piecewise"};
        double d;

        for (int i = 0; i < epsilonArray.size(); i++) {
            // double epsilon = epsilonArray[i];
            cout << "epsilon is: " << epsilonArray[i] << endl;
            for (int j = 0; j < 1; j++) {

                //cout << "weight number is: " << weightArray[j] << endl;
                for (int k = 0; k < penaltytypeArray.size(); k++) {
                    //create directories

                    //std::string directoryName = ".\\randomdist_it300_lam542_sym_epsilon_0.1_penaltytype_" + penaltytypeArray[k] + "_absolute_0.0to0.5\\";
                    //std::string directoryName = "..\\Calculations\\Clipped Random Initial Structure\\" + std::to_string(j) + "ClipRandomness_0.4to0.6_it400_lam542_sym_filterOff_periodicFalse_beta0_epsilon_0.1_penaltytype_" + penaltytypeArray[k] + "_0.0to0.5\\";
                    std::string directoryName = "E:\\Calculations\\Random Initial Structure\\periodicitydebug_it300_lam542_sym_filter2to3_periodicFalse_beta0_epsilon_0.5_penalty_piecewise0.0to0.5\\";
                    cout << "Storing data in : " << directoryName << endl;
                    std::filesystem::create_directories(directoryName);
                    std::filesystem::create_directories(directoryName + "/CoreStructure");
                    std::filesystem::create_directories(directoryName + "/E-Field");
                    std::filesystem::create_directories(directoryName + "/Model_output");
                    std::filesystem::create_directories(directoryName + "/Shape");
                    std::filesystem::create_directories(directoryName + "/ShapeSolid");

                    
                    VectorXi inputGeo;
                    VectorXd inputDiel;
                    int Nx, Ny, Nz;
                    int N = 0;
                    string pathCommonData = reader2.Get("Geometry", "pathCommonData", "UNKNOWN");
                    string pathPara = reader2.Get("Geometry", "pathPara", "UNKNOWN");
                    cout << "path common data file is: " << pathCommonData << endl;
                    tie(Nx, Ny, Nz, N, inputGeo, inputDiel) = getInputs(pathCommonData, pathPara);

                    // next 20 lines are if you want a random, symmetric distribution of pixel values

                    /*inputDiel = VectorXd::Zero(3 * N);
                    std::default_random_engine rnd{ std::random_device{}() };
                    std::uniform_real_distribution<double> dist(0.4, 0.6);

                    VectorXd inputdielxy = VectorXd::Zero(Nx * Ny);
                    int numxypixels = Nx * Ny;
                    for (int i = 0; i < numxypixels / 4; i++) {
                        int column = floor(i / (Nx / 2)) * Nx;
                        int row = i % (Nx / 2);
                        double r = dist(rnd);
                        //int index = i + (column * Nx / 2);
                        inputdielxy(column + row) = r;                   // original pixel
                        inputdielxy(column + Nx - row - 1) = r;       // reflection over y
                        inputdielxy(numxypixels - column - row - 1) = r;       // reflection over x and y
                        inputdielxy(numxypixels - Nx - column + row) = r;     // reflection over x
                    }

                    for (int i = 0; i < N; i++) {
                        double pixelvalue = inputdielxy(i % numxypixels);
                        inputDiel(3 * i) = pixelvalue;
                        inputDiel(3 * i + 1) = pixelvalue;
                        inputDiel(3 * i + 2) = pixelvalue;
                        
                    } */

                    string save_position = reader2.Get("Output", "saveDir", "UNKNOWN");       //output file
                    
                   // VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
                    /* for (int i = 0; i < total_space.size(); i++) {
                         cout << total_space(i) << endl;
                     } */


                    d = stod(reader2.Get("Grid", "d", "UNKNOWN"));

                    double lam = reader2.GetInteger("Input field", "lam", -1);
                    Vector3d n_K;
                    n_K(0) = reader2.GetFloat("Input field", "nKx", 0.0);
                    n_K(1) = reader2.GetFloat("Input field", "nKy", 0.0);
                    n_K(2) = reader2.GetFloat("Input field", "nKz", 0.0);
                    double E0 = 1.0;
                    Vector3d n_E0;
                    n_E0(0) = reader2.GetFloat("Input field", "nEx", 0.0);
                    n_E0(1) = reader2.GetFloat("Input field", "nEy", 0.0);
                    n_E0(2) = reader2.GetFloat("Input field", "nEz", 0.0);

                    list<string> mat_l{ reader2.Get("Material", "material1", "UNKNOWN"), reader2.Get("Material", "material2", "UNKNOWN") };
                    VectorXcd material = Get_X_material(mat_l, lam, "nm");

                    int MAX_ITERATION_DDA = reader2.GetInteger("DDA iteration", "MAX_ITERATION_DDA", -1);
                    int MAX_ITERATION_EVO = reader2.GetInteger("Evo Option", "MAX_ITERATION_EVO", -1);
                    double MAX_ERROR = reader2.GetFloat("DDA iteration", "MAX_ERROR", -1);

                    MatrixXi scope = find_scope_3_dim(&inputGeo);
                    int bindz = scope(2, 1) - scope(2, 0) + 1;
                    Vector3i bind(1, 1, bindz);

                    filterReader readTheFilter(reader2);
                    bool filter = readTheFilter.getFilter();
                    vector<filterinfo> filterList = readTheFilter.getFilterList();
                    int m, n;
                    int Lm, Ln;
                    m = reader2.GetInteger("Periodicity Option", "m", -1);
                    n = reader2.GetInteger("Periodicity Option", "n", -1);
                    Lm = reader2.GetInteger("Periodicity Option", "Lx", -1);
                    Ln = reader2.GetInteger("Periodicity Option", "Ly", -1);
                    FilterOption filterOpt(readTheFilter.getBetaMin(), readTheFilter.getBetaMax(), readTheFilter.getIta(), readTheFilter.getBetaType(), filterList);

                    symReader readSymmetry(reader2);
                    string symmetry = readSymmetry.getSymmetry();
                    vector<double> symAxis = readSymmetry.getSymAxis();
                    cout << "Periodicity is: " << readTheFilter.getPeriodic() << endl;
                    //StructureSpacePara structurespacepara(bind, &s, &inputDiel, filter, &filterOpt, symmetry, symAxis, readTheFilter.getPeriodic(), Lm, Ln); // line 1088 in SpacePara
                    StructureSpacePara structurespacepara(bind, &inputGeo, Nx, Ny, Nz, N, &inputDiel, filter, &filterOpt, symmetry, symAxis, readTheFilter.getPeriodic(), Lm, Ln);
                    cout << "SpacePara created" << endl;
                    CoreStructure CStr(&structurespacepara, d);
                    double nback = sqrt(real(material(0)));
                    cout << "CoreStructure created" << endl;
                    AProductCore Core(&CStr, lam, material, nback, m, n, Lm * d, Ln * d, "FCD");
                    cout << "AProductCore created" << endl;
                    DDAModel TestModel(&Core, n_K, E0, n_E0);
                    cout << "TestModel created" << endl;
                    ObjReader objReader(reader2);
                    string objName = objReader.GetObjName();
                    vector<double> objPara = objReader.GetObjPara();  //Focal spot position.

                    // CHANGING THIS TO TRUE TO SEE WHAT HAPPENS!!! IT WAS ORIGINALLY FALSE!!
                    bool HavePathRecord = false;
                    bool HaveOriginHeritage = false;
                    bool HaveAdjointHeritage = false;
                    double epsilonOld = reader2.GetFloat("Evo Option", "epsilon", 1.0);
                    cout << "epsilonOld is" << epsilonOld << "\n";



                    EvoDDAModel evoModel(objName, objPara, epsilonArray[i], HavePathRecord, HaveOriginHeritage, HaveAdjointHeritage, directoryName, &CStr, &TestModel);
                    evoModel.EvoOptimizationQuick(weightArray[j], penaltytypeArray[k], MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam"); // line 393 in EvoDDAModel
                }
            }
        }
        /*TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();

        string savePosition = reader2.Get("Output", "saveDir", "UNKNOWN");
        TestModel.output_to_file(savePosition + "Model_output\\", 0);

        string nameCommonData = save_position + "commondata.txt";
        ofstream Common;
        Common.open(save_position + "commondata.txt");
        Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
        Common << (spacepara.get_geometry()) << endl;
        Common << d << endl;
        Common << n_E0 << endl;
        Common << n_K << endl; */
        return;


    }
}

int main() {

    task();

}








