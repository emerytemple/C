#include "convert.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

double search_unit(struct unit* units, char *name)
{
	int i;
	double retval;

	for(i = 0; i < 4; ++i) {
		printf("%d, %s, %s\n", i, name, units[i].name);
		if(0 == strcmp(units[i].name, name)) {
			retval = units[i].amount;
			break;
		}
	}

	return retval;
}

bool isunittype(struct unit *units, char *name)
{
	int i;

	for(i = 0; i < 4; ++i) {
		if(0 == strcmp(units[i].name, name)) {
			return true;
		}
	}

	return false;
}

double Convert(double number, char* from_unit, char* to_unit)
{
	double from_val, to_val;
	struct unit *units;

	struct unit mass_units[] = {
		{"g",1000.0},
		{"kg",1.0},
		{"sl",0.0685218},
		{"lbm",2.20462}
	};

	struct unit distance_units[] = {
		{"mm",1000.0},
		{"cm",100.0},
		{"m",1.0},
		{"km",0.001},
		{"in",38.3701},
		{"ft",3.28084},
		{"yd",1.09361},
		{"mi",0.000621371}
	};

	struct unit time_units[] = {
		{"s",1.0},
		{"sec",1.0},
		{"min",0.0166667},
		{"hr",0.000277778}
	};

	struct unit pressure_units[] = {
		{"pa",1.0},
		{"atm",9.86923e-6},
		{"psi",0.000145038}
	};

	struct unit force_units[] = {
		{"n",1.0},
		{"lbf",0.224809}
	};

	struct unit angle_units[] = {
		{"deg",1.0},
		{"rad",M_PI/180.0}
	};

/*
	struct unit temperature_units[] = {
		{"c",},
		{"f",},
		{"k",},
		{"r",}
	};
*/

	if(true == isunittype(mass_units, from_unit)) {
		units = mass_units;
	}
	else if(true == isunittype(time_units, from_unit)) {
		units = time_units;
	}
	else if(true == isunittype(distance_units, from_unit)) {
		units = distance_units;
	}
	else if(true == isunittype(pressure_units, from_unit)) {
		units = pressure_units;
	}
	else if(true == isunittype(force_units, from_unit)) {
		units = force_units;
	}
	else if(true == isunittype(angle_units, from_unit)) {
		units = angle_units;
	}
	else {
	}

	from_val = search_unit(units, from_unit);
	to_val = search_unit(units, to_unit);

	return number*(to_val/from_val);
}
