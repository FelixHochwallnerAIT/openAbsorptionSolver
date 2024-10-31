// Defining custom functions

// Include necessary modules
#include <vector>
#include <algorithm>

// Define necessary functions

// Find next lower value of point in data
scalar findLower(std::vector<scalar> data, scalar point)
{
	scalar lower=data[0];
	for(const auto& value : data)
	{
		if(value<=point)
		{
			lower = value;
		}
	}
	return lower;
}

// Find next higher value of point in data
scalar findHigher(std::vector<scalar> data, scalar point)
{
	scalar higher=data.back();
	for(const auto& value : data)
	{
		if(value>=point)
		{
			higher = value;
			break;
		}
	}
	return higher;
}

// Find index of point in data
int findIndex(std::vector<scalar> data, scalar point)
{
	int idx = -1;
	std::vector<scalar>::iterator it = find(data.begin(), data.end(), point);
	if (it != data.end()) // Point found
	{
		idx = distance(data.begin(),it);
	}
	else // Point not found
	{
		idx = -1;
	}
	return idx;
}

// Linear interpolate between x1 and x2
scalar linearInterpolation(scalar x1, scalar f_x1, scalar x2, scalar f_x2, scalar x)
{
	scalar result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
	return result;
}
