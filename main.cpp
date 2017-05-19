#include <ExPDE.hpp>
#include "math.h"
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <Datatypes/Types.hpp>

#include <Datatypes/Coordinate.hpp>
#include <Datatypes/Variable/variable.hpp>
#include <IO/VTK.hpp>
#include <Datatypes/Variable/variable_typedefs.hpp>


// boundary condition

class MarkerFunctor
{
public:
	using value_type = bool;

	MarkerFunctor( const double Ly_, const double hy_ )
	    : Ly( Ly_ )
	    , hy( hy_ )

	{
	}

	value_type evaluate( double y )
	{
		if ( y  <  hy )
            return true;
		return false;
	}

private:
	double Ly;
	double hy;

};


//call in the directory of tsk2.exe:
//tsk2 lambda ppw cfl N_timesteps

int main(int argc , char *argv[])
{
	// namespace of the ExPDE library
	using namespace ExPDE;

	//source parameters
    //assume there is no loss, namely sigma = 0
    //the wave propagates in vacuum

	double epsl = 8.854187817e-12; // vacuum permittivity
	double miu = ( M_PI * 4.0e-7 ); // magnetic constant
	double c = 299792458.0; // light speed in vacuum
	double lambda = 0.550e-6; //wavelength
	double omega = 2 * M_PI * c / lambda;

	//grid parameters

	const double Lx = 10 * lambda; //domain size
	const double Ly = 5 * lambda;
	double ppw = 10;
	const double hx = lambda / ppw; //mesh size in x and y direction
	const double hy = lambda / ppw;
    const double Nx = Lx / hx;//number of grid points in x direction
	const double Ny = Ly / hy;//number of grid points in y direction
	double cfl = 1;
	double timestep =  cfl / ( c * sqrt(1. / (hx * hx) + 1. / (hy * hy) ) ) ;
	int N_timestep = 100; //iteration times

	//call and input parameters

	switch ( argc )
	{
		case 5:
			lambda     = std::atof(argv[1]);
			ppw        = std::atof(argv[2]);
			cfl        = std::atof(argv[3]);
			N_timestep = std::atoi(argv[4]);
			break;

		default:
			std::cerr << "programm-call should be:"
			          << "\"./[ProgsName] lambda ppw cfl timesteps\"!" << std::endl;
			return EXIT_FAILURE;
	}

	//create grid

	Staggered_grid grid( Nx, Ny, 1, {0., 0., 0.}, {Lx, Ly, hx}, periodic, periodic_no, periodic_no );

	// boundary condition functor

	MarkerFunctor mfunc( Ly, hy );
	ET::Functor< MarkerFunctor, bool, double> marker_functor( mfunc );

    //create a Coordinate iterators for the staggered grid to be used with the functor

	Coordinate< _X_ > X( grid );
	Coordinate< _Y_ > Y( grid );

	//create staggered E and not-staggered H, initialized with all 0

	Variable< double, not_staggered, staggered, not_staggered > Hx( grid, 0.0 );
	Variable< double, staggered, not_staggered, not_staggered > Hy( grid, 0.0 );
	Variable< double, not_staggered, not_staggered, not_staggered >Ez( grid, 0.0 );
	Variable< double, not_staggered, not_staggered, not_staggered >epsl_( grid, epsl );
	Variable< bool, not_staggered, staggered, not_staggered >bound( grid, marker_functor(Y) );

	//iteration

    for (int i = 0; i < N_timestep; ++i)
	{
        double t = i * timestep; //physical time

        //update three components

        Hx = bound | 1.0 * sin(omega * t);
        Ez = Ez + (timestep / 2) / epsl * ( ( E(Hy) - W(Hy) ) / hx - ( N(Hx) - S(Hx) ) / hy );
        Hy = Hy + timestep / miu * ( E(Ez) - W(Ez) ) / hx;
        Hx = Hx - timestep / miu * ( N(Ez) - S(Ez) ) / hy;

        IO::VTK::XMLPrinter( Hx, IO::concat_string( "Hx_", i ) );
		IO::VTK::XMLPrinter( Hy, IO::concat_string( "Hy_", i ) );
		IO::VTK::XMLPrinter( Ez, IO::concat_string( "Ez_", i ) );
	}

	return 0;
}
