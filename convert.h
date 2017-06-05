#ifndef _CONVERT_H_
#define _CONVERT_H_

struct unit {
	char *name;
	double amount;
};

double Convert(double number, char *from_unit, char *to_unit);

#endif
