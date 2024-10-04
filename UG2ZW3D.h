#pragma once

#ifndef UG2ZW3D_H
#define UG2ZW3D_H
#endif

/*
本头文件由邹锋编写
联系方式QQ609719845（微信同号）    E-mail:609719845@qq.com
周末有空会时更新
可能程序中会存在错误BUG，可以指正一下
如果你也写有其他的可以发我
本源码目的是一个人忙碌，换取更多人无缝转换代码。节约更多时间
一人写一个，500人就写500个，更好
*/

#include <windows.h>
#include <fstream>
#include <list>
#include <shlobj.h>
#include <tchar.h>
#include <atlstr.h>	  
#include <strstream>
#include <io.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <time.h> 
 #include <map>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <stdlib.h>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm> 
#include <cmath> //C++
#include <math.h> //C语言
//fabs()计算绝对值
//sqrt() 计算平方根
//cbrt()计算立方根
//pow()幂运算
//ceil() 向上取整
//floor() 向下取整
using namespace std;
#include "api/inc/VxApi.h"

static double PI = 3.1415926535897932384626433832795;
 typedef unsigned int tag_t;
#define UF_MODL_CYLINDRICAL_FACE        16   /* UF_MODL_ask_face_type */
#define UF_MODL_CONICAL_FACE            17   /* UF_MODL_ask_face_type */
#define UF_MODL_SPHERICAL_FACE          18   /* UF_MODL_ask_face_type */
#define UF_MODL_TOROIDAL_FACE           19   /* UF_MODL_ask_face_type */
#define UF_MODL_SWEPT_FACE              20   /* UF_MODL_ask_face_type */
#define UF_MODL_PLANAR_FACE             22   /* UF_MODL_ask_face_type */
#define UF_MODL_BLENDING_FACE           23   /* UF_MODL_ask_face_type */
#define UF_MODL_PARAMETRIC_FACE         43   /* UF_MODL_ask_face_type */
#define UF_MODL_OFFSET_FACE             65   /* UF_MODL_ask_face_type */
#define UF_MODL_FOREIGN_FACE            66   /* UF_MODL_ask_face_type */

#define UF_MODL_LINEAR_EDGE             3001    /* UF_MODL_ask_edge_type */
#define UF_MODL_CIRCULAR_EDGE           3002    /* UF_MODL_ask_edge_type */
#define UF_MODL_ELLIPTICAL_EDGE         3003    /* UF_MODL_ask_edge_type */
#define UF_MODL_INTERSECTION_EDGE       3004    /* UF_MODL_ask_edge_type */
#define UF_MODL_SPLINE_EDGE             3005    /* UF_MODL_ask_edge_type */

#define UF_MODL_OPEN_CURVE               0  /* UF_MODL_ask_curve_periodicity*/
#define UF_MODL_CLOSED_PERIODIC_CURVE    1  /* UF_MODL_ask_curve_periodicity*/



 enum UF_MODL_boolean_body
 {
	 UF_MODL_TARGET_BODY = 0,
	 UF_MODL_TOOL_BODY = 1
 };
 typedef enum UF_MODL_boolean_body UF_MODL_boolean_body_e_t;

 enum UF_FEATURE_SIGNS
 {
	 UF_NULLSIGN = 0,    /* create new target solid */
	 UF_POSITIVE = 1,    /* add to target solid */
	 UF_NEGATIVE = 2,     /* subtract from target solid */
	 UF_UNSIGNED = 3,    /* intersect with target solid */

/* the following are new settings introduced for the function
   UF_MODL_ask_feature_boolean
   and are intended to eventually replace the previous settings. */

	UF_NO_BOOLEAN = 4,  /* feature has not been booleaned */
	UF_TOP_TARGET = 5,  /* feature is the "top target" feature, it has no
						  "parent" features but does have tool features */
	UF_UNITE = 6,       /* feature has been united to target solid */
	UF_SUBTRACT = 7,    /* feature has been subtracted from target solid */
	UF_INTERSECT = 8,   /* feature has been intersected with target solid */
	UF_DEFORM_POSITIVE = 9, /* feature used to deform the positive side
							of the target sheet */
	UF_DEFORM_NEGATIVE = 10 /* feature used to deform the negative side
							of the target sheet */
 };
 typedef enum UF_FEATURE_SIGNS UF_FEATURE_SIGN;


 struct UF_CURVE_line_s {
	 double      start_point[3];  /* line start point */
	 double      end_point[3];    /* line end point */
 };

 typedef struct UF_CURVE_line_s UF_CURVE_line_t,
	 * UF_CURVE_line_p_t;


 struct UF_CURVE_arc_s {
	 int       matrix_tag;             /* matrix for the CSYS the arc is in */
	 double      start_angle;            /* expressed in radians */
	 double      end_angle;              /* expressed in radians */
	 double      arc_center[3];       /* center of the arc */
	 double      radius;                 /* radius of the arc */
 };

 typedef struct UF_CURVE_arc_s UF_CURVE_arc_t,
	 * UF_CURVE_arc_p_t;
	 



 //转换点OK
 static void ZwMathMapPoint(svxMatrix FormMat, svxPoint FormPt, svxMatrix ToMat, svxPoint& ToPt)
 {
	 FormMat.identity = 0;//表明初始化的矩阵已经被修改了
	 ToMat.identity = 0;//表明初始化的矩阵已经被修改了

	 svxPoint pt = FormPt;//相对于矩阵B的坐标
	 svxMatrix invertFormMat;
	 cvxMatInvert(&FormMat, &invertFormMat);//求相对坐标系矩阵A的逆
	 svxMatrix transMat;
	 cvxMatMult(&invertFormMat, &ToMat, &transMat);//矩阵相乘
	 cvxPntTransform(&transMat, &pt);//点转换
	 ToPt = pt;
	 //princ(pt);
 }
	 
//句柄列表转ID数组
static vector <int> hand2id(int count, szwEntityHandle* handleList)
{
	vector <int> EntIds;
	int* pEntIds = nullptr;
	ZwMemoryAlloc(count * sizeof(int), (void**)&pEntIds);
	ZwEntityIdGet(count, handleList, pEntIds);
	for (size_t i = 0; i < count; i++)
	{
		EntIds.push_back(pEntIds[i]);
		//princ(pEntIds[i]);
	}
	//ZwMemoryFree((void**)pEntIds);
	cvxMemFree((void**)&pEntIds);

	return EntIds;
}

//单个句柄转ID
static int hand2id(szwEntityHandle handleid)
{
	szwEntityHandle handleList[1] = { handleid };
	int pEntIds[1] = { 0 };
	ZwEntityIdGet(1, handleList, pEntIds);
	return pEntIds[0];
}

//句柄要释放 ZwEntityHandleFree(&entityHandle);
//单个句柄转ID
static int id2hand(int objid, szwEntityHandle* entityHandle)
{
	int ret = ZwEntityIdTransfer(1, &objid, entityHandle);
	return ret;
}


//移除参数
static int zwdelparam(int objectid)
{
	int objid = 0;
	int ret = cvxPartDefeature(1, &objectid);
	return objid;
}


//找体的面
static int pk_ask_body_faces(szwEntityHandle shapeHandle, vector <szwEntityHandle>& bodyfaces)
{
	bodyfaces.clear();
	int count;
	szwEntityHandle* faceList;
	ZwShapeFaceListGet(shapeHandle, &count, &faceList);
	for (size_t i = 0; i < count; i++)
	{
		bodyfaces.push_back(faceList[i]);
	}
	ZwEntityHandleListFree(count, &faceList);
	return 0;
}


//PK_BODY_unite_bodies(pk_cyl1, 1, booer_body, &n_bodies, &bodies);
//求和1和多
static	int PK_BODY_unite_bodies(int pk_target_body, int num_tool, int* pk_tool_body, int* n_bodies, int** bodies)
{
	cvxPartBool(VX_BOOL_ADD, pk_target_body, num_tool, pk_tool_body, 0);
	return 0;
}

//求和 1和1
static	int PK_Unite_Bodies(int pk_target_body, int pk_tool_body)
{
	cvxPartBool(VX_BOOL_ADD, pk_target_body, 1, &pk_tool_body, 0);
	return 0;
}

//求差  1和多
static	int PK_BODY_subtract_bodies(int pk_target_body, int num_tool, int* pk_tool_body, int* n_bodies, int** bodies)
{
	cvxPartBool(VX_BOOL_REMOVE, pk_target_body, num_tool, pk_tool_body, 0);
	return 0;
}

//求差 1和1
static	int PK_Subtract_Bodies(int pk_target_body, int pk_tool_body)
{
	cvxPartBool(VX_BOOL_REMOVE, pk_target_body, 1, &pk_tool_body, 0);
	return 0;
}


static	int UF_MODL_ask_face_type(int face_id, int* facetype)
{
	svxSrfPrim SrfPrim;
	cvxPartInqFaceSrfPrim(face_id, &SrfPrim);
	if (SrfPrim.srfType == VX_SF_PRIM_PLN)
	{
		*facetype = UF_MODL_PLANAR_FACE;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_SPH)
	{
		*facetype = UF_MODL_SPHERICAL_FACE;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_CON)
	{
		*facetype = UF_MODL_CONICAL_FACE;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_CYL)
	{
		*facetype = UF_MODL_CYLINDRICAL_FACE;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_NURB)
	{
		*facetype = UF_MODL_PARAMETRIC_FACE;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_ELLSO)
	{
		*facetype = UF_MODL_FOREIGN_FACE;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_TORUS)
	{
		//圆环面
		*facetype = UF_MODL_TOROIDAL_FACE;
	}
	cvxSrfPrimFree(&SrfPrim);
	return 0;
}

static	int UF_MODL_ask_edge_type(int edge_id, int* edgetype)
{
	svxCurve Crv;
	cvxPartInqCurve(edge_id, 1, &Crv);
	if (VX_CRV_LINE == Crv.Type)
	{
		*edgetype = UF_MODL_LINEAR_EDGE;
	}
	else if (VX_CRV_ARC == Crv.Type || VX_CRV_CIRCLE == Crv.Type)
	{
		*edgetype = UF_MODL_CIRCULAR_EDGE;
	}
	else if (VX_CRV_NURB == Crv.Type)
	{
		*edgetype = UF_MODL_SPLINE_EDGE;
	}
	cvxCurveFree(&Crv);
	return 0;
}

//注意跟UG与ZW是相反的
static	int UF_MODL_ask_curve_periodicity(int curve_id, int* isopen)
{
	//0: the specified edge is not an open edge
	//- 1 : the specified edge is an open edge
	int ret = cvxPartInqEdgeOpen(curve_id);
	if (ret)
	{
		*isopen = 0;
	}
	else
	{
		*isopen = 1;
	}
	return 0;
}


static	int UF_MODL_ask_list_count(vector <int> objects, int* num)
{
	*num = (int)objects.size();
	return 0;
}

static	int UF_MODL_ask_list_item(vector <int> objects, int index, int* objid)
{
	if (index < objects.size())
	{
		*objid = objects[index];
	}
	else
	{
		return 1;
	}
	return 0;
}

static	int UF_MODL_put_list_item(vector <int>& objects, int objid)
{
	objects.push_back(objid);
	return 0;
}

static	int UF_MODL_ask_face_body(int face, int* bodyid)
{
	szwEntityHandle facehand;
	id2hand(face, &facehand);
	szwEntityHandle shape;
	*bodyid = hand2id(shape);
	//cvxPartInqFaceShape(face, bodyid);
	return 0;
}

static	int UF_MODL_ask_body_faces(int bodyid, vector <int>& objects)
{
	objects.clear();
	int count;
	int* faces;
	cvxPartInqShapeFaces(bodyid, &count, &faces);
	for (size_t j = 0; j < count; j++)
	{
		objects.push_back(faces[j]);
	}
	cvxMemFree((void**)&faces);
	return 0;
}

static	int UF_MODL_ask_body_edges(int bodyid, vector <int>& objects)
{
	objects.clear();
	int eCount;
	int* Edegs = NULL;
	cvxPartInqShapeEdges(bodyid, &eCount, &Edegs);
	for (size_t c = 0; c < eCount; c++)
	{
		objects.push_back(Edegs[c]);
	}
	cvxMemFree((void**)&Edegs);
	return 0;
}

static	int UF_MODL_ask_face_edges(int face_tag, vector <int>& faceedges)
{
	faceedges.clear();
	int eCount;
	int* Edegs = NULL;
	cvxPartInqFaceEdges(face_tag, &eCount, &Edegs);
	for (size_t c = 0; c < eCount; c++)
	{
		faceedges.push_back(Edegs[c]);
	}
	cvxMemFree((void**)&Edegs);
	return 0;
}

static	int UF_MTX3_initialize(const double x_vec[3], const double y_vec[3], double mtx[9])
{
	svxPoint p0 = { 0 };
	svxVector v1 = { x_vec[0],x_vec[1],x_vec[2] };
	svxVector v2 = { y_vec[0],y_vec[1],y_vec[2] };
	svxMatrix matx;
	cvxMatPntVecs(&p0, &v1, &v2, &matx);
	mtx[0] = matx.xx;
	mtx[1] = matx.xy;
	mtx[2] = matx.xz;
	mtx[3] = matx.yx;
	mtx[4] = matx.yy;
	mtx[5] = matx.yz;
	mtx[6] = matx.zx;
	mtx[7] = matx.zy;
	mtx[8] = matx.zz;
	return 0;
}

static	int UF_MTX3_initialize_x(const double x_vec[3], const double y_vec[3], double(&mtx)[9])
{
	svxPoint p0 = { 0 };
	svxVector v1 = { x_vec[0],x_vec[1],x_vec[2] };
	svxMatrix matx;
	cvxMatPntVec(&p0, &v1, &matx);
	mtx[6] = matx.xx * -1;
	mtx[7] = matx.xy * -1;
	mtx[8] = matx.xz * -1;
	mtx[3] = matx.yx;
	mtx[4] = matx.yy;
	mtx[5] = matx.yz;
	mtx[0] = matx.zx;
	mtx[1] = matx.zy;
	mtx[2] = matx.zz;
}

static	int UF_MTX3_initialize_z(const double x_vec[3], const double y_vec[3], double(&mtx)[9])
{
	svxPoint p0 = { 0 };
	svxVector v1 = { x_vec[0],x_vec[1],x_vec[2] };
	svxMatrix matx;
	cvxMatPntVec(&p0, &v1, &matx);
	mtx[0] = matx.xx;
	mtx[1] = matx.xy;
	mtx[2] = matx.xz;
	mtx[3] = matx.yx;
	mtx[4] = matx.yy;
	mtx[5] = matx.yz;
	mtx[6] = matx.zx;
	mtx[7] = matx.zy;
	mtx[8] = matx.zz;
	return 0;
}

static	int UF_VEC3_copy(const double vec[3], double(&copy_vec)[3])
{
	for (size_t i = 0; i < 3; i++)
	{
		copy_vec[i] = vec[i];
	}
	return 0;
}

static void UF_VEC3_cross(const double vec1[3], const double ve2[3], double(&cross_vec)[3])
{
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	szwVector v2 = { ve2[0],ve2[1],ve2[2] };
	szwVector v3;
	ZwVectorCrossProduct(v1, v2, &v3);
	cross_vec[0] = v3.x;
	cross_vec[1] = v3.y;
	cross_vec[2] = v3.z;
}

static void UF_VEC3_dot(const double vec1[3], const double ve2[3], double* dotPorduct)
{
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	szwVector v2 = { ve2[0],ve2[1],ve2[2] };
	ZwVectorDot(v1, v2, dotPorduct);
}

//计算某个点按指定某个point2的方向移动dist距离后的点pnt_on_seg
static void UF_VEC3_convex_comb(double dist, const double start_pt[3], const double point2[3], double(&pnt_on_seg)[3])
{
	svxPoint p1 = { start_pt[0],start_pt[1],start_pt[2] };
	svxPoint p2 = { point2[0],point2[1],point2[2] };
	svxVector v1;
	cvxVecInit(&p1, &p2, &v1);
	cvxPntTranslate(&p1, &v1, dist);
	pnt_on_seg[0] = p1.x;
	pnt_on_seg[1] = p1.y;
	pnt_on_seg[2] = p1.z;

}

static void UF_VEC3_angle_between(const double vec1[3], const double ve2[3], double vec_ccw[3], double* angle)
{
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	szwVector v2 = { ve2[0],ve2[1],ve2[2] };
	ZwVectorAngleGet(v1, v2, angle);
}

static void UF_VEC3_unitize(const double vec1[3], double tolerance, double* magnitude, double unit_vec[3])
{
	*magnitude = 0;
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	ZwVectorNormalize(&v1);
	unit_vec[0] = v1.x;
	unit_vec[1] = v1.y;
	unit_vec[2] = v1.z;
}

static void  UF_VEC3_ask_perpendicular(const double vec1[3], double(&perpendicular_vec)[3])
{
	szwVector v0 = { vec1[0],vec1[1],vec1[2] };
	szwVector v1 = { 0 };
	ZwVectorPerpendicularGet(v0, &v1);
	perpendicular_vec[0] = v1.x;
	perpendicular_vec[1] = v1.y;
	perpendicular_vec[2] = v1.z;
}

static void UF_VEC3_is_perpendicular(const double vec1[3], const double vec2[3], const double tol, int* is_perp)
{
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	szwVector v2 = { vec2[0],vec2[1],vec2[2] };
	ZwVectorIsPerpendicular(v1, v2, tol, is_perp);

}

static void UF_VEC3_is_equal(const double vec1[3], const double vec2[3], const double tol, int* isEqual)
{
	szwVector vector1 = { vec1[0],vec1[1],vec1[2] };
	szwVector vector2 = { vec2[0],vec2[1],vec2[2] };
	ZwVectorIsEqual(vector1, vector2, tol, isEqual);
}

static void UF_VEC3_is_parallel(const double vec1[3], const double vec2[3], const double tol, int* ifParallel)
{
	szwVector vector1 = { vec1[0],vec1[1],vec1[2] };
	szwVector vector2 = { vec2[0],vec2[1],vec2[2] };
	ZwVectorParallelCheck(0, vector1, vector2, tol, ifParallel);
}

static void UF_VEC3_sub(const double vec1[3], const double vec2[3], double(&vec_diff)[3])
{
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	szwVector v2 = { vec2[0],vec2[1],vec2[2] };
	szwVector v3 = { 0 };
	ZwVectorSubtraction(v1, v2, &v3);
	vec_diff[0] = v3.x;
	vec_diff[1] = v3.y;
	vec_diff[2] = v3.z;
}

static void UF_VEC3_add(const double vec1[3], const double vec2[3], double(&vec_sum)[3])
{
	szwVector v1 = { vec1[0],vec1[1],vec1[2] };
	szwVector v2 = { vec2[0],vec2[1],vec2[2] };
	szwVector v3 = { 0 };
	ZwVectorSum(v1, v2, &v3);
	vec_sum[0] = v3.x;
	vec_sum[1] = v3.y;
	vec_sum[2] = v3.z;
}

static void UF_VEC3_triple(const double vec1[3], const double vec2[3], const double vec3[3], double* tripleScaleProduct)
{
	szwVector vector1 = { vec1[0],vec1[1],vec1[2] };
	szwVector vector2 = { vec2[0],vec2[1],vec2[2] };
	szwVector vector3 = { vec3[0],vec3[1],vec3[2] };
	ZwVectorTripleScaleProduct(vector1, vector2, vector3, tripleScaleProduct);
}

static void UF_VEC3_mag(const double vec1[3], double* vectorMagnitude)
{
	szwVector vector1 = { vec1[0],vec1[1],vec1[2] };
	ZwVectorMagnitude(vector1, vectorMagnitude);
}

static void UF_VEC3_scale(double scale, const double vec1[3], double(&scaled_vec)[3])
{
	szwVector vector1 = { vec1[0],vec1[1],vec1[2] };
	szwVector scaleVector;
	ZwVectorScaling(scale, vector1, &scaleVector);
	scaled_vec[0] = scaleVector.x;
	scaled_vec[1] = scaleVector.y;
	scaled_vec[2] = scaleVector.z;
}

static void UF_VEC3_midpt(const double pt1[3], const double pt2[3], double(&cenpt)[3])
{
	for (size_t i = 0; i < 3; i++)
	{
		cenpt[i] = (pt1[i] + pt2[i]) / 2;
	}
}

//测距离
static int UF_VEC3_distance(double pt1[3], double pt2[3], double* distance)
{
	//svxPoint P1 = { pt1[0] ,pt1[1],pt1[2] };
	//svxPoint P2 = { pt2[0] ,pt2[1],pt2[2] };
	//*distance = cvxPntDist(&P1, &P2);
	double x1 = pt1[0];
	double x2 = pt2[0];
	double y1 = pt1[1];
	double y2 = pt2[1];
	double z1 = pt1[2];
	double z2 = pt2[2];
	//\[D = \sqrt{ (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2 } \]
	*distance = std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2) + std::pow(z2 - z1, 2));
	return 0;
}

static int UF_MTX3_x_vec(double matx[9], double(&dir)[3])
{
	for (int k = 0; k < 3; k++)
	{
		dir[k] = matx[k];
	}
	return 0;
}
static int UF_MTX3_y_vec(double matx[9], double(&dir)[3])
{
	for (int k = 0; k < 3; k++)
	{
		dir[k] = matx[k + 3];
	}
	return 0;
}
static int UF_MTX3_z_vec(double matx[9], double(&dir)[3])
{
	for (int k = 0; k < 3; k++)
	{
		dir[k] = matx[k + 6];
	}
	return 0;
}
static int UF_VEC3_negate(double dir[3], double(&dir1)[3])
{
	for (int k = 0; k < 3; k++)
	{
		dir1[k] = -1 * dir[k];
	}
	return 0;
}

//移除参数
static void UF_MODL_delete_object_parm(int objectid)
{
	zwdelparam(objectid);

}

//移除参数
static void UF_MODL_delete_object_parm1(int objectid)
{
	zwdelparam(objectid);
}

static int UF_MODL_create_cyl1(UF_FEATURE_SIGN sign, double origin[3], char* height, char* diam, double direction[3], int* cyl_obj_id)
{
	svxPoint pt = { origin[0], origin[1], origin[2] };
	svxVector dir = { direction[0], direction[1], direction[2] };
	svxCylData Cyl;
	cvxPartCylInit(&Cyl);
	Cyl.axis.Pnt = pt;
	Cyl.axis.Dir = dir;
	Cyl.useAxis = 1;
	Cyl.Radius = atof(diam);
	Cyl.Length = atof(height);
	Cyl.Center = pt;

	//UF_NULLSIGN = 0,    /* create new target solid */
	//UF_POSITIVE = 1,    /* add to target solid */
	//UF_NEGATIVE = 2,     /* subtract from target solid */
	//UF_UNSIGNED = 3,    /* intersect with target solid */

	//VX_BOOL_NONE = 0, /**< @brief base */
	//VX_BOOL_ADD = 1, /**< @brief add */
	//VX_BOOL_REMOVE = 2, /**< @brief remove */
	//VX_BOOL_INTERSECT = 3  /**< @brief intersect */

	Cyl.Combine = (evxBoolType)sign;
	int ret = cvxPartCyl(&Cyl, cyl_obj_id);
	return ret;
}

static int UF_MODL_create_cylinder(UF_FEATURE_SIGN sign, tag_t targ_tag, double origin[3], char* height, char* diam, double direction[3], int* cyl_obj_id)
{
	UF_MODL_create_cyl1(UF_NULLSIGN, origin, height, diam, direction, cyl_obj_id);
	int ret = cvxPartBool((evxBoolType)sign, targ_tag, 1, cyl_obj_id, 0);
	return ret;
}

//另外自己写了一个
static int UF_MODL_create_block2(double point[3], double mtx[9], double box_size[3], int* blk_obj_id)
{
	svxMatrix matx;
	svxPoint Cenpt = { point[0],point[1],point[2] };
	svxVector xvec = { mtx[0],mtx[1],mtx[2] };
	svxVector yvec = { mtx[3],mtx[4],mtx[5] };
	svxVector zvec = { 0,0,0 };

	cvxMatPntVecs(&Cenpt, &xvec, &yvec, &matx);
	cvxMatGetPntVecs(&matx, &Cenpt, &xvec, &yvec, &zvec);

	int idPlane;

	svxPlaneData Plane;
	Plane.method = VX_PLANE_DYNAMIC;
	Plane.inpUnion.dynamic.PosPnt = Cenpt;
	Plane.inpUnion.dynamic.XAxis = xvec;
	Plane.inpUnion.dynamic.YAxis = yvec;
	Plane.dOffset = 0;
	int ret = cvxPartPlaneNew(&Plane, &idPlane);

	cvxPntTranslate(&Cenpt, &xvec, box_size[0] / 2);
	cvxPntTranslate(&Cenpt, &yvec, box_size[1] / 2);
	cvxPntTranslate(&Cenpt, &zvec, box_size[2] / 2);

	svxBoxData Box;
	cvxPartBoxInit(&Box);
	cvxMemZero((void*)&Box, sizeof(Box));
	strcpy_s(Box.ftrName, sizeof(Box.ftrName), "");
	Box.useAxis = 0;
	Box.idPlane = idPlane;
	Box.Center = Cenpt;
	Box.X = box_size[0];
	Box.Y = box_size[1];
	Box.Z = box_size[2];
	Box.Combine = VX_BOOL_NONE;
	int iRet = cvxPartBox(&Box, blk_obj_id);
	return iRet;
}

//注意生成方块它是取中心点，坐标系可以 输入，要自己修改成自己想要的
static int UF_MODL_create_block1(UF_FEATURE_SIGN sign, double origin[3], char* edge_len[3], int* blk_obj_id)
{
	//ZwLCSMatrixGet(szwMatrix * matrix);
	//要获取当前的坐标系

	svxPoint pt = { origin[0], origin[1], origin[2] };
	svxBoxData blkbox;
	cvxPartBoxInit(&blkbox);
	blkbox.Center = pt;
	blkbox.X = atof(edge_len[0]);
	blkbox.Y = atof(edge_len[1]);
	blkbox.Z = atof(edge_len[2]);
	int ret = cvxPartBox(&blkbox, blk_obj_id);
	return ret;
}

static int UF_MODL_create_block(UF_FEATURE_SIGN sign, tag_t targ_tag, double origin[3], char* edge_len[3], int* blk_obj_id)
{
	UF_MODL_create_block1(UF_NULLSIGN, origin, edge_len, blk_obj_id);
	int ret = cvxPartBool((evxBoolType)sign, targ_tag, 1, blk_obj_id, 0);
	return ret;
}

static int UF_MODL_create_cone1(UF_FEATURE_SIGN sign, double origin[3], char* height, char* diam[2], double direction[3], int* cone_obj_id)
{
	svxPoint pt = { origin[0], origin[1], origin[2] };
	svxVector dir = { direction[0], direction[1], direction[2] };
	svxConeData Con;
	cvxPartConeInit(&Con);
	Con.axis.Pnt = pt;
	Con.axis.Dir = dir;
	Con.useAxis = 1;
	Con.Radius1 = atof(diam[0]);
	Con.Radius2 = atof(diam[1]);
	Con.Length = atof(height);
	Con.Center = pt;
	int ret = cvxPartCone(&Con, cone_obj_id);
	return ret;
}

static int UF_MODL_create_cone(UF_FEATURE_SIGN sign, tag_t targ_tag, double origin[3], char* height, char* diam[2], double direction[3], int* cone_tag)
{
	UF_MODL_create_cone1(UF_NULLSIGN, origin, height, diam, direction, cone_tag);
	int ret = cvxPartBool((evxBoolType)sign, targ_tag, 1, cone_tag, 0);
	return ret;
}

static int UF_MODL_create_sphere1(UF_FEATURE_SIGN sign, double origin[3], char* diam, int* sphere_obj_id)
{
	svxPoint pt = { origin[0], origin[1], origin[2] };
	svxSphereData Sphere;
	cvxPartSphereInit(&Sphere);
	Sphere.Center = pt;
	Sphere.Combine = (evxBoolType)sign;
	Sphere.Radius = atof(diam);
	int ret = cvxPartSphere(&Sphere, sphere_obj_id);
	return ret;
}

static int UF_MODL_create_sphere(UF_FEATURE_SIGN sign, tag_t targ_tag, double origin[3], char* diam, int* sphere_tag)
{
	UF_MODL_create_sphere1(UF_NULLSIGN, origin, diam, sphere_tag);
	int ret = cvxPartBool((evxBoolType)sign, targ_tag, 1, sphere_tag, 0);
	return ret;
}

static int UF_CURVE_create_point(double pt[3], int* pt_tag)
{
	svxPoint p1 = { pt[0], pt[1], pt[2] };
	int ret = cvxPartPnt(&p1, pt_tag);
	return ret;
}

static int UF_CURVE_create_line(UF_CURVE_line_t line_coords, int* line_tag)
{
	svxPoint p1 = { line_coords.start_point[0],line_coords.start_point[1],line_coords.start_point[2] };
	svxPoint p2 = { line_coords.end_point[0],line_coords.end_point[1],line_coords.end_point[2] };
	int ret = cvxPartLine2pt(&p1, &p2, line_tag);
	return ret;
}

static int UF_CURVE_ask_line_data(int line_tag, double(&start_point)[3], double(&end_point)[3])
{
	szwEntityHandle entityHandle;
	id2hand(line_tag, &entityHandle);
	int ret = 0;
	int isCurve;
	ZwEntityCurveCheck(entityHandle, &isCurve);
	if (isCurve)
	{
		szwPoint start, end;
		ret = ZwCurveEndPointGet(entityHandle, &start, &end);
		start_point[0] = start.x;
		start_point[1] = start.y;
		start_point[2] = start.z;

		end_point[0] = end.x;
		end_point[1] = end.y;
		end_point[2] = end.z;
	}
	ZwEntityHandleFree(&entityHandle);
	return ret;
}

//UF_CURVE_line_t line_coords;
//UF_CURVE_ask_line_data(line_tag, &line_coords);
static int UF_CURVE_ask_line_data(int line_tag, UF_CURVE_line_t* line_coords)
{
	szwEntityHandle entityHandle;
	id2hand(line_tag, &entityHandle);
	szwPoint start, end;
	int ret = ZwCurveEndPointGet(entityHandle, &start, &end);
	line_coords->start_point[0] = start.x;
	line_coords->start_point[1] = start.y;
	line_coords->start_point[2] = start.z;

	line_coords->end_point[0] = end.x;
	line_coords->end_point[1] = end.y;
	line_coords->end_point[2] = end.z;

	ZwEntityHandleFree(&entityHandle);
	return ret;
}

//建议用3点画圆比较好
static int UF_CURVE_create_arc(UF_CURVE_arc_t arc_coords, int* arc_tag)
{
	int ret = 0;
	if (fabs(arc_coords.end_angle - arc_coords.start_angle - PI) < 0.001)
	{
		svxCircleData arcData;
		cvxPartCircleInit(&arcData);
		arcData.type = VX_CIRCLE_RADIUS;
		arcData.idAlignPln = arc_coords.matrix_tag;//最好使用平面矩阵，在平面中心
		arcData.useDiameter = 0;// (0: radius; 1: diameter)
		arcData.radOrDia = arc_coords.radius;
		ret = cvxPartCircle(&arcData, arc_tag);
		//cvxPartCir3pt()//建议用3点画圆比较好
	}
	else
	{
		svxArcData arcData;
		cvxPartArcInit(&arcData);
		arcData.arcType = VX_ARC_ANGLE;
		arcData.idAlignPln = arc_coords.matrix_tag;//最好使用平面矩阵，在平面中心
		arcData.radius = arc_coords.radius;
		arcData.startAngle = arc_coords.start_angle;
		arcData.arcAngle = arc_coords.end_angle - arc_coords.start_angle;//ZW可能使用起始角度与圆弧角度，不是弧度，也不是终止角度
		ret = cvxPartArc(&arcData, arc_tag);
		//cvxPartArc3pt();//建议用3点画圆弧比较好
	}
	return ret;
}

//获取曲线基中一节的长度，也可以 总长
static int UF_CURVE_ask_arc_length(int curve, double start_param, double end_param, int unit_flag, double* length)
{
	szwEntityHandle curveHandle;
	id2hand(curve, &curveHandle);
	int ret = ZwCurveSegmentLengthGet(curveHandle, start_param, end_param, length);
	//cvxCrvLen2(int idCurve, double T1, double T2, double *Length);
	ZwEntityHandleFree(&curveHandle);
	return ret;
}

//获取曲线总长度
static int UF_CURVE_ask_arc_length(int curve, double* length)
{
	szwEntityHandle curveHandle;
	id2hand(curve, &curveHandle);
	int ret = ZwCurveLengthGet(curveHandle, length);
	ZwEntityHandleFree(&curveHandle);
	return ret;
}


//获取曲线总长度
static int UF_MODL_ask_curve_props(int curve, double parm, double(&point)[3], double(&tangent)[3], double p_norm[3], double b_norm[3], double* torsion, double* rad_of_cur)
{
	szwEntityHandle curveHandle;
	id2hand(curve, &curveHandle);
	szwPoint p1;
	//int ret = ZwCurvePointGetByLengthFraction(curveHandle, parm,&p1);//这个只获取到点
	szwCurveDerivative evaluate;
	szwVector normal;
	int ret = ZwCurveDifferentiate(curveHandle, parm, 3, &evaluate, &normal);

	point[0] = p1.x;
	point[1] = p1.y;
	point[2] = p1.z;

	tangent[0] = normal.x;
	tangent[1] = normal.y;
	tangent[2] = normal.z;
	ZwEntityHandleFree(&curveHandle);
}

//获取曲线质心
static int UF_CURVE_ask_centroid(int curve, double(&point)[3])
{
	szwEntityHandle curveHandle;
	id2hand(curve, &curveHandle);
	szwPoint centroidPoint;
	int ret = ZwCurveCentroidPointGet(curveHandle, &centroidPoint);

	point[0] = centroidPoint.x;
	point[1] = centroidPoint.y;
	point[2] = centroidPoint.z;
	ZwEntityHandleFree(&curveHandle);
	return ret;
}


/*
int curve;
int ret = cvxGetEnt("选择对象", evxEntInpOpt::VX_INP_CURVE, 1, &curve);
int numpts;
double* pts=NULL;
UF_MODL_ask_curve_points(curve, 0.01, 0.01, 2.5, &numpts, &pts);
int ss = 0;
for (size_t i = 0; i < numpts*3;)
{
	double pt[3] = { 0 };
	pt[0] = pts[0 + i];
	pt[1] = pts[1 + i];
	pt[2] = pts[2 + i];
	int pttag;
	svxPoint p1 = { pt[0], pt[1], pt[2] };
	cvxPartPnt(&p1, &pttag);
	i = i + 3;
	//char msg[256] = "";
	//sprintf_s(msg, "%f %f %f \n", pt[0], pt[1], pt[2]);
	//cvxMsgDisp(msg);
}
free(pts);
*/
//获取曲线上的点
//要释放点free(pts);
static int UF_MODL_ask_curve_points(int curve, double ctol, double atol, double steptol, int* n_pts, double** pts)
{
	szwEntityHandle curveHandle;
	id2hand(curve, &curveHandle);
	szwPoint* points;
	int ret = ZwCurveTessellationPointListGet(curveHandle, ctol, steptol, n_pts, &points);
	if (*n_pts > 0)
	{
		*pts = (double*)malloc(sizeof(double) * (*n_pts * 3));
		if (*pts == NULL) {
			cvxMsgDisp("UF_MODL_ask_curve_points Memory allocation failed");
			return 1;
		}
		int n = 0;
		for (size_t i = 0; i < *n_pts; i++)
		{
			(*pts)[n] = points[i].x;
			(*pts)[n + 1] = points[i].y;
			(*pts)[n + 2] = points[i].z;
			n = n + 3;
		}
	}
	ZwMemoryFree((void**)&points);
	ZwEntityHandleFree(&curveHandle);
	return ret;
}

//使用vector
static int UF_MODL_ask_curve_points(int curve, double ctol, double atol, double steptol, int* n_pts, vector<double> pts)
{
	szwEntityHandle curveHandle;
	id2hand(curve, &curveHandle);
	szwPoint* points;
	int ret = ZwCurveTessellationPointListGet(curveHandle, ctol, steptol, n_pts, &points);
	if (*n_pts > 0)
	{
		for (size_t i = 0; i < *n_pts; i++)
		{
			pts.push_back(points[i].x);
			pts.push_back(points[i].y);
			pts.push_back(points[i].z);
		}
	}
	ZwMemoryFree((void**)&points);
	ZwEntityHandleFree(&curveHandle);
	return ret;
}

static int UF_OBJ_delete_object(int obj)
{
	szwEntityHandle entityHandle;
	id2hand(obj, &entityHandle);
	int ret = ZwEntityDelete(entityHandle);
	ZwEntityHandleFree(&entityHandle);
	return ret;
}

static int UF_OBJ_delete_array_of_objects(const int num_objects, int* objectid, int** status)
{
	*status = 0;
	szwEntityHandle* entityHandles;
	entityHandles = (szwEntityHandle*)malloc(sizeof(szwEntityHandle) * num_objects);
	ZwEntityIdTransfer(num_objects, objectid, entityHandles);
	int ret = ZwEntityListDelete(num_objects, entityHandles);
	for (int i = 0; i < num_objects; ++i)
	{
		ZwEntityHandleFree(&entityHandles[i]);
	}
	free(entityHandles);
	return ret;
}


static int UF_OBJ_ask_color(int obj, int* colorvalue)
{

}

static int UF_OBJ_set_color(int obj, int colorvalue)
{
	//UG216种颜色的RGB值
	int UGColor[217][3] =
	{ 0,0,0
	,255,255,255
	,255,250,191
	,255,245,167
	,255,242,131
	,255,240,126
	,255,255,0
	,225,245,255
	,215,255,199
	,204,255,153
	,204,255,102
	,204,255,51
	,255,232,98
	,189,255,255
	,153,255,202
	,242,242,242
	,221,235,214
	,255,240,169
	,220,227,154
	,153,255,255
	,114,255,205
	,220,239,237
	,208,232,213
	,193,193,144
	,207,222,107
	,78,255,255
	,192,227,223
	,175,219,219
	,203,230,200
	,51,255,51
	,188,237,52
	,0,255,255
	,195,224,241
	,159,213,210
	,180,219,190
	,153,255,102
	,0,255,0
	,249,215,231
	,255,215,210
	,255,204,153
	,251,199,160
	,255,192,76
	,255,202,0
	,202,199,228
	,204,204,204
	,240,204,134
	,255,220,58
	,204,204,51
	,230,219,73
	,207,215,232
	,223,223,223
	,232,236,216
	,227,211,157
	,225,194,72
	,158,202,61
	,102,204,255
	,96,205,210
	,195,225,190
	,227,232,134
	,169,227,0
	,102,204,0
	,153,204,255
	,0,207,207
	,171,214,160
	,156,208,126
	,160,212,168
	,136,187,112
	,152,201,235
	,154,207,168
	,0,205,134
	,183,207,195
	,110,179,97
	,120,192,140
	,255,169,255
	,245,212,185
	,255,153,153
	,246,160,73
	,255,153,51
	,255,153,0
	,179,195,224
	,221,218,201
	,207,154,156
	,201,151,101
	,204,153,51
	,171,152,82
	,154,186,214
	,151,149,197
	,153,153,153
	,196,191,165
	,155,144,0
	,178,177,103
	,137,182,255
	,102,153,204
	,132,173,155
	,55,170,136
	,141,148,52
	,118,159,42
	,51,153,255
	,0,157,217
	,37,159,165
	,0,150,111
	,51,153,51
	,72,144,106
	,0,153,255
	,0,176,240
	,0,153,153
	,0,172,101
	,94,133,63
	,0,153,0
	,245,191,216
	,255,179,180
	,255,144,144
	,255,102,102
	,246,160,104
	,255,102,0
	,181,181,255
	,211,159,200
	,204,102,153
	,247,175,132
	,230,123,18
	,210,107,55
	,192,210,225
	,173,168,212
	,177,177,177
	,153,102,102
	,153,102,51
	,111,104,0
	,114,106,255
	,102,102,204
	,152,170,175
	,102,102,102
	,102,102,51
	,78,97,0
	,130,182,196
	,51,102,204
	,54,96,146
	,0,88,40
	,51,102,51
	,74,115,27
	,65,96,255
	,49,133,155
	,0,102,153
	,0,102,102
	,28,73,3
	,0,102,0
	,255,109,255
	,255,58,225
	,255,48,48
	,255,139,139
	,217,86,76
	,234,110,165
	,211,189,171
	,204,51,204
	,182,146,110
	,179,129,93
	,204,51,51
	,156,116,41
	,209,179,158
	,166,163,196
	,125,125,125
	,153,51,102
	,153,51,51
	,102,51,0
	,198,168,141
	,102,51,204
	,118,106,150
	,81,31,127
	,109,45,10
	,73,69,46
	,51,51,255
	,152,176,216
	,128,162,180
	,62,51,84
	,51,51,51
	,41,36,34
	,0,68,255
	,56,120,192
	,32,88,103
	,0,63,94
	,0,71,64
	,0,58,0
	,255,0,255
	,255,81,81
	,177,117,147
	,144,50,106
	,238,149,189
	,255,0,0
	,134,115,97
	,192,48,164
	,164,102,0
	,144,80,7
	,147,139,100
	,172,55,19
	,122,115,181
	,92,73,128
	,95,80,69
	,65,55,50
	,126,0,0
	,153,0,0
	,102,0,255
	,102,0,204
	,76,76,76
	,66,25,66
	,101,0,68
	,102,0,0
	,86,177,255
	,0,57,154
	,42,73,114
	,57,17,85
	,36,0,36
	,12,12,12
	,0,0,255
	,0,0,192
	,36,64,97
	,18,40,109
	,15,36,62
	,0,0,0
	};

	svxColor rgb;
	rgb.r = UGColor[colorvalue][0];
	rgb.g = UGColor[colorvalue][1];
	rgb.b = UGColor[colorvalue][2];
	int ret = cvxEntRgbSet(rgb, 1, &obj);
	return ret;
}

static int UF_OBJ_ask_blank_status(int obj, int* blanked)
{
	szwEntityHandle entity;
	id2hand(obj, &entity);
	int ret = ZwEntityBlankGet(entity, blanked);
	ZwEntityHandleFree(&entity);
	return ret;
}

static int UF_OBJ_set_blank_status(int obj, int blanked)
{
	szwEntityHandle entity;
	id2hand(obj, &entity);
	int ret = ZwEntityBlankSet(1, &entity, blanked);
	ZwEntityHandleFree(&entity);
	return ret;
}

static int UF_OBJ_ask_name(int obj, char* name)
{
	szwEntityHandle entity;
	id2hand(obj, &entity);
	char nameTag[256] = "";
	int size = 256;
	int ret = ZwEntityNameGet(entity, size, nameTag);
	ZwEntityHandleFree(&entity);
	strcpy_s(name,256, nameTag);
	return ret;
}

static int UF_OBJ_set_name(int obj, char* name)
{
	szwEntityHandle entity;
	id2hand(obj, &entity);
	int ret = ZwEntityNameSet(entity, name);
	ZwEntityHandleFree(&entity);
	return ret;
}

