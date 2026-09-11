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

#ifndef __pointGeom_hpp__
#define __pointGeom_hpp__

#include "types.hpp"
#include "dictionary.hpp"
#include "iIstream.hpp"
#include "iOstream.hpp"

namespace pFlow
{

class pointGeom
{
private:

	realx3 	point_{0,0,0};

public:

	TypeInfoNV("point");

		INLINE_FUNCTION_HD
		pointGeom() = default;

		INLINE_FUNCTION_HD
		explicit pointGeom(const realx3& p)
		:
			point_(p)
		{}

		FUNCTION_H
		explicit pointGeom(const dictionary& dict);

		FUNCTION_H
		explicit pointGeom(iIstream& is);

		INLINE_FUNCTION_HD
		pointGeom(const pointGeom&) = default;

		INLINE_FUNCTION_HD
		pointGeom(pointGeom&&) = default;

		INLINE_FUNCTION_HD
		pointGeom& operator=(const pointGeom&) = default;

		INLINE_FUNCTION_HD
		pointGeom& operator=(pointGeom&&) = default;

		INLINE_FUNCTION_HD
		~pointGeom()=default;

		INLINE_FUNCTION_HD
		bool isInside(const realx3& p)const
		{
			return equal(p, point_);
		}

		INLINE_FUNCTION_HD
		const realx3& minPoint()const
		{
			return point_;
		}

		INLINE_FUNCTION_HD
		const realx3& maxPoint()const
		{
			return point_;
		}

		INLINE_FUNCTION_HD
		const realx3& center()const
		{
			return point_;
		}

		INLINE_FUNCTION_HD
		real volume()const
		{
			return static_cast<real>(0);
		}

		FUNCTION_H
		bool read(iIstream & is);

		FUNCTION_H
		bool write(iOstream& os)const;

		FUNCTION_H
		bool read(const dictionary& dict);

		FUNCTION_H
		bool write(dictionary& dict)const;
};

FUNCTION_H
iIstream& operator >>(iIstream& is, pointGeom& p);

FUNCTION_H
iOstream& operator << (iOstream& os, const pointGeom& p);

INLINE_FUNCTION_HD
bool equal(const pointGeom& p1, const pointGeom& p2, real tol = smallValue)
{
	return equal(p1.center(), p2.center(), tol);
}

INLINE_FUNCTION_HD
bool operator ==(const pointGeom& p1, const pointGeom& p2)
{
	return equal(p1, p2);
}

}

#endif
