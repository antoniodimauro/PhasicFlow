/*------------------------------- PhasicFlow ---------------------------------
  Copyright (C): Antonio Di Mauro
  email: antoniodimauro03@gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of PhasicFlow, a CFD-DEM stack built on phasicFlow and
  PhasicFlowPlus (www.cemf.ir). It is free software: you can redistribute it
  and/or modify it under the terms of the GNU General Public License v3 or
  any later version.

  It is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.
-----------------------------------------------------------------------------*/

#include "pointGeom.hpp"

FUNCTION_H
pFlow::pointGeom::pointGeom
(
	const dictionary & dict
)
:
	point_
	(
		dict.getVal<realx3>("point")
	)
{}

FUNCTION_H
pFlow::pointGeom::pointGeom
(
	iIstream& is
)
{
	if( !read(is))
	{
		ioErrorInFile(is.name(), is.lineNumber())<<
		"error in reading point from file. \n";
		fatalExit;
	}
}

FUNCTION_H
bool pFlow::pointGeom::read(iIstream & is)
{
	if(!is.nextData<realx3>("point", point_)) return false;
	return true;
}

FUNCTION_H
bool pFlow::pointGeom::write(iOstream& os)const
{
	os.writeWordEntry("point", point_);
	return os.check(FUNCTION_NAME);
}

FUNCTION_H
bool pFlow::pointGeom::read
(
	const dictionary& dict
)
{
	point_ = dict.getVal<realx3>("point");
	return true;
}

FUNCTION_H
bool pFlow::pointGeom::write
(
	dictionary& dict
)const
{
	if(!dict.add("point", point_))
	{
		fatalErrorInFunction<<
		"  error in writing point to dictionary "<<dict.globalName()<<endl;
		return false;
	}

	return true;
}

FUNCTION_H
pFlow::iIstream& pFlow::operator >>(iIstream& is, pointGeom& p)
{
	if(! p.read(is))
	{
		ioErrorInFile(is.name(), is.lineNumber())<<
		"error in reading point. \n";
		fatalExit;
	}
	return is;
}

FUNCTION_H
pFlow::iOstream& pFlow::operator <<(iOstream& os, const pointGeom& p)
{

	if(! p.write(os))
	{
		ioErrorInFile(os.name(), os.lineNumber())<<
		"error in writing point. \n";
		fatalExit;
	}
	return os;
}
