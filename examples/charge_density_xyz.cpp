#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdarg>
#include <sstream>
#include <cctype>

using namespace std;

// ---------------------------- //
// Functions
// ---------------------------- //

// FUNCTION TO SAVE THE FILE INTO AN ARRAY
void read_from_the_file(vector<string>& vector_file,string file_name);

// EXTRACT THE VARIABLES Nx, Ny, Nz, number_atoms AND CREATE A 1ROW-VECTOR WITH ALL THE CHARGE DENSITY
void extract_variables(vector<string>& vector_file,vector<float>& charge_density, int& nx, int& ny, int& nz, int& number_atoms);

// CREATE COLOR FOR CHARGE DENSITY
string create_color(float charge);

// CREATE AN OUPTUT XYZ FILE
void create_output_file(vector<float> charge_density,string file_name,int nx, int ny, int nz);

// ---------------------------- //
// Main
// ---------------------------- //

int main(int argc, char* argv[])
{

    string file_name;
    string part_nb;
    int number_points; // number of lines - 2
    vector<string> vector_file;
    vector<float> charge_density;
    
    int nx(0),ny(0),nz(0); // number of points in each direction of the grid
    int number_atoms(0);
    
    cout << "Name of the file to read (with extension) : " << endl;
    cin >> file_name;
    
    read_from_the_file(vector_file,file_name);
    
    number_points = vector_file.size()-2;
    
    cout << "Number of points : " << number_points << endl;
    
    extract_variables(vector_file,charge_density,nx,ny,nz,number_atoms);
    
    create_output_file(charge_density,file_name,nx,ny,nz);
}

// ---------------------------- //
// Definition of functions
// ---------------------------- //

void read_from_the_file(vector<string>& vector_file,string file_name)
{
    std::ifstream file(file_name.c_str());
    if(file.is_open())
    {
        string new_line;
        while (getline(file,new_line))
        {
            vector_file.push_back(new_line);
        }
    }
    else cout << "Unable to open file";
    file.close();
}

void extract_variables(vector<string>& vector_file,vector<float>& charge_density, int& nx, int& ny, int& nz, int& number_atoms)
{
    std::string str("");
    std::string str0("");
    
    const char *cstr("") ;
    
    float d1(0.0),d2(0.0),d3(0.0),d4(0.0),d5(0.0);
    
    // read the second line and attribute the values to nx,ny,nz and number_atoms
    str = vector_file[1];
    cstr = str.c_str();
    sscanf(cstr,"%d %d %d %d %d %d %d",&nx,&ny,&nz,&nx,&ny,&nz,&number_atoms);

    int start(0);
    start = 5+number_atoms;
    
    int incr(0);
    // create the vector of densities
    for (unsigned int i(start); i<vector_file.size(); i++)
    {
        
        str = vector_file[i];
        cstr = str.c_str();
        
        sscanf(cstr, "%E %E %E %E %E",&d1,&d2,&d3,&d4,&d5);
        charge_density.push_back(d1);
        incr = incr+1;
        charge_density.push_back(d2);
        incr = incr+1;
        charge_density.push_back(d3);
        incr = incr+1;
        charge_density.push_back(d4);
        incr = incr+1;
        charge_density.push_back(d5);
        incr = incr+1;
    
    }

}

string create_color(float charge)
{
    if (charge <= 0.0005)
    {
        return "a";
    }
    
    else if((charge > 0.0005) && (charge <= 0.001))
    {
        return "b";
    }
    
    else if((charge > 0.001) && (charge <=0.01))
    {
        return "c";
    }
    
    else if((charge > 0.01) && (charge <= 0.05))
    {
        return "d";
    }
    
    else if((charge > 0.05) && (charge <= 0.1))
    {
        return "e";
    }
    
    else if((charge > 0.1) && (charge <= 0.5))
    {
        return "f";
    }
    else if((charge > 0.5))
    {
        return "g";
    }

    else
        return 0;
}

void create_output_file(vector<float> charge_density,string file_name,int nx, int ny, int nz)
{
    
    string output_name;
    output_name = file_name + "_.xyz";
    std::ofstream file_out;

    file_out.open(output_name.c_str());
    
    file_out << charge_density.size() << endl;
    file_out << "colored " << endl;
    
    int incr(0);
    
    for (unsigned int k(1); k<=nz; k++)
    {
        for (unsigned int j(1); j<=ny; j++)
        {
            for (unsigned int i(1); i<=nx; i++)
            {
                file_out << create_color(charge_density[incr]) << " " << i << ". " << j << ". " << k << ". " << charge_density[incr] << endl;
                incr = i + (j-1) * nx + (k-1)*nx*ny;
            }
        }
    }
    
    
    file_out.close();
    
}









