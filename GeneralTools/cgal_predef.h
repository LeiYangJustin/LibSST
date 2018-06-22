#ifndef CGAL_PREDEF_H
#define CGAL_PREDEF_H

// delaunay triangulation
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>  DelaunayT;
typedef K::Point_2 KPoint2;
typedef K::Vector_2 KVector2;


#include <CGAL/Polygon_2.h>
typedef CGAL::Polygon_2<K>						KPolygon2;

#include <CGAL/Segment_2.h>
typedef CGAL::Segment_2<K>						KSegment2;


#endif // ! CGAL_PREDEF_H

