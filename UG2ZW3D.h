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


#define _USE_MATH_DEFINES
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
#include <iostream>
#include <cmath>
//fabs()计算绝对值
//sqrt() 计算平方根
//cbrt()计算立方根
//pow()幂运算
//ceil() 向上取整
//floor() 向下取整
using namespace std;

//#include "uf_object_types.h"
#include "D:\\ZWAPI\\uf_object_types.h"
#include "api/inc/VxApi.h"

//static double PI = 3.1415926535897932384626433832795;
#ifndef TAG_T_DEFINED
#define TAG_T_DEFINED
typedef unsigned int tag_t;
typedef tag_t* tag_p_t;
#endif

#ifndef UF_LIST_T_DEFINED
#define UF_LIST_T_DEFINED
typedef struct uf_list_s* uf_list_p_t;
struct uf_list_s {
	tag_t                  eid;  /* Object ID */
	struct uf_list_s* next;  /* Pointer to the next OID in the list */
};

typedef struct uf_list_s uf_list_t;
#endif

#ifndef NULL_TAG
#define NULL_TAG      ((tag_t)0)
#endif

#ifndef LOGICAL_DEFINED
#define LOGICAL_DEFINED
#if defined(__cplusplus)
typedef bool logical;
#else
typedef unsigned char logical;
#endif
#endif

#ifndef BYTE_DEFINED
#define BYTE_DEFINED
typedef unsigned char byte;
#endif

#if !defined(true) && !defined(__cplusplus)
#define true   1
#define false  0
#endif

#ifndef TRUE
#define TRUE   1
#define FALSE  0
#endif


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
	 
#define UF_ATTR_MAX_STRING_LEN 132
#define UF_ATTR_MAX_TITLE_LEN   50

 /*****************************************************************************
 * Attribute type definitions
 ****************************************************************************/
#define UF_ATTR_integer    1
#define UF_ATTR_real       2
#define UF_ATTR_time       3
#define UF_ATTR_null       4
#define UF_ATTR_string     5
#define UF_ATTR_any        6
#define UF_ATTR_reference  7

 /*****************************************************************************
 * Attribute value
 ****************************************************************************/
 union UF_ATTR_value_u
 {
	 int    integer;    /* If the attribute is an integer attribute,
						   this can be used to access the value. */
	 double real;       /* If the attribute is a floating point attribute,
						   this can be used to access the value. */
	 int    time[2];    /* If the attribute is a date/time attribute,
						   this can be used to access the value.
						   time[0] contains the NX computational date
						   time[1] contains the NX computational time. */
	 char* string;    /* If the attribute is a string attribute, then
						   this is a pointer to the value. The maximum length
						   of this string is limited by UF_ATTR_MAX_STRING_LEN */

	 char* reference; /* If the attribute is a reference attribute, then
						   this is a pointer to the value.  The value may
						   have an embedded expression with the following
						   syntax:
							   <Xm.n@exp_name> or <Xm,n@exp_name>
						   The X indicates that an expression is being
						   referenced.  The m specifies the minimum field
						   width.  If necessary, it is padded on the left to
						   make up the field width.  The n specifies the
						   maximum number of digits after the decimal point of
						   the expression value.  The decimal point can be
						   specified as either . or , thus either m.n or m,n
						   are acceptable.  */
 };

 typedef union UF_ATTR_value_u UF_ATTR_value_u_t, * UF_ATTR_value_u_p_t;
 /*****************************************************************************
 * Typed attribute value
 ****************************************************************************/
 struct UF_ATTR_value_s
 {
	 int type;                 /* The type of the attribute.  Valid types are:
								  UF_ATTR_integer
								  UF_ATTR_real
								  UF_ATTR_time
								  UF_ATTR_null
								  UF_ATTR_string
								  UF_ATTR_reference
							   */
	 UF_ATTR_value_u_t value;  /* The attribute value */
 };

 typedef struct UF_ATTR_value_s UF_ATTR_value_t, * UF_ATTR_value_p_t;


 static void princ(svxMatrix matx)
 {
	 char msg[256] = "";
	 sprintf_s(msg, "matx pt= %f,  %f,  %f", matx.xt, matx.yt, matx.zt);
	 cvxMsgDisp(msg);
	 sprintf_s(msg, "matx X= %f,  %f,  %f", matx.xx, matx.xy, matx.xz);
	 cvxMsgDisp(msg);
	 sprintf_s(msg, "matx Y= %f,  %f,  %f", matx.yx, matx.yy, matx.yz);
	 cvxMsgDisp(msg);
	 sprintf_s(msg, "matx Z= %f,  %f,  %f", matx.zx, matx.zy, matx.zz);
	 cvxMsgDisp(msg);
 }

 static void princ_matx(double matx[9])
 {
	 char msg[1024] = "";
	 sprintf_s(msg, "matx X= %f, %f, %f", matx[0], matx[1], matx[2]);
	 cvxMsgDisp(msg);
	 sprintf_s(msg, "matx Y= %f, %f, %f", matx[3], matx[4], matx[5]);
	 cvxMsgDisp(msg);
	 sprintf_s(msg, "matx Z= %f, %f, %f", matx[6], matx[7], matx[8]);
	 cvxMsgDisp(msg);
 }

 static void princ(string msg)
 {
	 char msg1[1024];
	 strcpy_s(msg1, msg.c_str());
	 cvxMsgDisp(msg1);
 }


 static void princ(char msg[256])
 {
	 cvxMsgDisp(msg);
 }

 static void princ(double value1)
 {
	 char msg[256] = "";
	 sprintf_s(msg, "%f", value1);
	 cvxMsgDisp(msg);
 }
 static void princ(double point1[3])
 {
	 char msg[1024] = "";
	 sprintf_s(msg, "%f %f %f", point1[0], point1[1], point1[2]);
	 cvxMsgDisp(msg);
 }

 static void princ(int value1)
 {
	 char msg[256] = "";
	 sprintf_s(msg, "%d", value1);
	 cvxMsgDisp(msg);
 }

 static void princ(char* canshu, int value1)
 {
	 char msg[256] = "";
	 sprintf_s(msg, canshu, value1);
	 cvxMsgDisp(msg);
 }


 static void princ(string canshu, string value1)
 {
	 char msg[256] = "";
	 sprintf_s(msg, canshu.c_str(), value1.c_str());
	 cvxMsgDisp(msg);
 }

 static void princ(size_t value1)
 {
	 char msg[256] = "";
	 sprintf_s(msg, "%d", (int)value1);
	 cvxMsgDisp(msg);
 }

 static void princ(svxVector vec)
 {
	 char msg[256] = "";
	 sprintf_s(msg, "%f %f %f", vec.x, vec.y, vec.z);
	 cvxMsgDisp(msg);
 }
 static void princ(char* canshu, svxVector vec)
 {
	 char msg[256] = "";
	 char canshu1[256] = "";
	 strcpy_s(canshu1, canshu);
	 strcat_s(canshu1, "%f, %f, %f");
	 sprintf_s(msg, canshu1, vec.x, vec.y, vec.z);
	 cvxMsgDisp(msg);
 }


 static void princ(svxPoint point)
 {
	 char msg[256] = "";
	 sprintf_s(msg, "%f %f %f", point.x, point.y, point.z);
	 cvxMsgDisp(msg);
 }
 static void princ(char* canshu, svxPoint point)
 {
	 char msg[256] = "";

	 char canshu1[256] = "";
	 strcpy_s(canshu1, canshu);
	 strcat_s(canshu1, "%f, %f, %f");
	 sprintf_s(msg, canshu1, point.x, point.y, point.z);
	 cvxMsgDisp(msg);
 }

 static void uc1601(char* msg, int option)
 {
	 if (option==0)
	 {
		 cvxMsgDisp(msg);
	 }
	 else if (option == 1)
	 {
		 cvxGetResponse(1, msg);
	 }	 
 }

static void UF_UI_open_listing_window()
 {
 }
 
 static void UF_UI_write_listing_window(const char* msg)
 {
	 cvxMsgDisp(msg);
 }

 // 创建一个链表
 static int UF_MODL_create_list(uf_list_p_t* list) {
	 *list = NULL;
	 return 0;
 }

 // 查询链表数量
 static int UF_MODL_ask_list_count(uf_list_p_t list, int* count) {
	 *count = 0;
	 uf_list_t* current = list;
	 while (current != NULL) {
		 (*count)++;
		 current = current->next;
	 }
	 return 0;
 }

 // 查询链表对象
 static int UF_MODL_ask_list_item(uf_list_p_t list, int index, tag_t* object) {
	 uf_list_t* current = list;
	 int i = 0;
	 while (current != NULL && i < index) {
		 current = current->next;
		 i++;
	 }
	 if (current != NULL) {
		 *object = current->eid;
	 }
	 else {
		 *object = 0; // 表示索引超出范围
	 }
	 return 0;
 }


 // 把对象添加到链表中
 static int UF_MODL_put_list_item(uf_list_p_t& list, tag_t obj_id) {
	 uf_list_t* newNode = (uf_list_t*)malloc(sizeof(uf_list_t));
	 if (newNode == NULL) {
		 fprintf(stderr, "Memory allocation failed\n");
		 return 1;
	 }
	 newNode->eid = obj_id;
	 newNode->next = NULL;

	 if (list == NULL) {
		 list = newNode;
	 }
	 else {
		 uf_list_t* current = list;
		 while (current->next != NULL) {
			 current = current->next;
		 }
		 current->next = newNode;
	 }
	 return 0;
 }


 // 从链表中删除对象
 static  int UF_MODL_delete_list_item(uf_list_s** head, tag_t data) {
	 uf_list_s* current = *head;
	 uf_list_s* previous = NULL;

	 // 如果要删除的是头节点
	 if (current != NULL && current->eid == data) {
		 *head = current->next;
		 free(current);
		 return 0; 
	 }

	 // 遍历链表找到要删除的节点
	 while (current != NULL && current->eid != data) {
		 previous = current;
		 current = current->next;
	 }

	 // 如果找到了要删除的节点
	 if (current != NULL) {
		 // 从链表中移除它
		 previous->next = current->next;
		 free(current);
	 }
	 return 0;
 }

 //删除链表并清空整个链表的函数
 static int UF_MODL_delete_list(uf_list_p_t* head) {
	 uf_list_p_t current = *head;
	 uf_list_p_t next;

	 while (current != NULL) {
		 next = current->next; // 保存下一个节点的指针
		 free(current);        // 释放当前节点的内存
		 current = next;       // 移动到下一个节点
	 }

	 *head = NULL; // 将头指针设为 NULL
	 return 0;
 }



 static void UF_MTX3_point(svxMatrix matx, svxPoint& point)
 {
	 point.x = matx.xt;
	 point.y = matx.yt;
	 point.z = matx.zt;
 }

 static void UF_MTX3_x_vec(svxMatrix matx, svxVector& vec)
 {
	 vec.x = matx.xx;
	 vec.y = matx.xy;
	 vec.z = matx.xz;
 }



 static void UF_MTX3_y_vec(svxMatrix matx, svxVector& vec)
 {
	 vec.x = matx.yx;
	 vec.y = matx.yy;
	 vec.z = matx.yz;
 }

 static void UF_MTX3_z_vec(svxMatrix matx, svxVector& vec)
 {
	 vec.x = matx.zx;
	 vec.y = matx.zy;
	 vec.z = matx.zz;
 }



 static	void UF_MTX3_x_vec(svxMatrix mtx, double xvec[3])
 {
	 xvec[0] = mtx.xx;
	 xvec[1] = mtx.xy;
	 xvec[2] = mtx.xz;
 }

 static	void UF_MTX3_y_vec(svxMatrix mtx, double xvec[3])
 {
	 xvec[0] = mtx.yx;
	 xvec[1] = mtx.yy;
	 xvec[2] = mtx.yz;
 }

 static	void UF_MTX3_z_vec(svxMatrix mtx, double zvec[3])
 {
	 zvec[0] = mtx.zx;
	 zvec[1] = mtx.zy;
	 zvec[2] = mtx.zz;
 }

 static	void UF_MTX3_point(svxMatrix mtx, double point[3])
 {
	 point[0] = mtx.xt;
	 point[1] = mtx.yt;
	 point[2] = mtx.zt;
 }


 static	void UF_VEC3_copy(svxPoint pt, double copy_vec[3])
 {
	 copy_vec[0] = pt.x;
	 copy_vec[1] = pt.y;
	 copy_vec[2] = pt.z;
 }

 static	void UF_VEC3_copy(double copy_vec[3], svxPoint& pt)
 {
	 pt.x = copy_vec[0];
	 pt.y = copy_vec[1];
	 pt.z = copy_vec[2];
 }


 static	void UF_VEC3_copy(svxVector vec, double(&copy_vec)[3])
 {
	 copy_vec[0] = vec.x;
	 copy_vec[1] = vec.y;
	 copy_vec[2] = vec.z;
 }


 static	void UF_VEC3_copy(double copy_vec[3], svxVector& vec)
 {
	 vec.x = copy_vec[0];
	 vec.y = copy_vec[1];
	 vec.z = copy_vec[2];
 }



 //绝对坐标转用户坐标
 //UF_CSYS_map_point(UF_CSYS_WORK_COORDS, abs_point,UF_CSYS_ROOT_WCS_COORDS, wcs_point);
 static void UF_CSYS_map_point(svxMatrix from_csys, svxPoint from_pt, svxMatrix to_csys, svxPoint& to_pt)
 {
	 svxMatrix to_csys_invert, transMat;
	 cvxMatInit(&transMat);
	 cvxMatInvert(&to_csys, &to_csys_invert);
	 cvxMatMult(&from_csys, &to_csys_invert, &transMat);
	 cvxPntTransform(&transMat, &from_pt);
	 to_pt = from_pt;
 }


 //参数1和3这个改成了对象ID
 static void UF_CSYS_map_point(int from_csys_id, svxPoint from_pt, int to_csys_id, svxPoint& to_pt)
 {
	 svxCSYSData from_csys_data, to_csys_data;
	 cvxCSYSGetData(from_csys_id, &from_csys_data);
	 cvxCSYSGetData(to_csys_id, &to_csys_data);
	 UF_CSYS_map_point(from_csys_data.Frame, from_pt, to_csys_data.Frame, to_pt);
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
static int zwdelparam(int& objectid)
{
	int* label = NULL;
	cvxEntLabelGet(objectid, &label);
	//移除参数
	int count = 1;//造型数量
	cvxPartDefeature(count, &objectid);

	//重新获取圆柱ID
	cvxEntByLabel(label, 1, &objectid);
	cvxMemFree((void**)&label);

	return objectid;
}
//移除参数
static void zwdelparam(int Count, int* objectid)
{
	int ret = cvxPartDefeature(Count, objectid);
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
	//szwEntityHandle facehand;
	//id2hand(face, &facehand);
	//szwEntityHandle shape;
	//*bodyid = hand2id(shape);
	cvxPartInqFaceShape(face, bodyid);
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


static int UF_VEC3_is_parallel1(const double a1[3], const double a2[3], double a3, int* result)
{
	double v7, v8, v9, v10, v11, v12;

	if (fabs(a1[0]) >= 1.0e19 || fabs(a1[1]) >= 1.0e19 || fabs(a1[2]) >= 1.0e19)
	{
		return 1;
	}
	if (fabs(a2[0]) >= 1.0e19 || fabs(a2[1]) >= 1.0e19 || fabs(a2[2]) >= 1.0e19)
	{
		return 1;
	}

	v7 = a1[0];
	*result = 1;
	if (fabs(a1[1] * a1[1] + v7 * v7 + a1[2] * a1[2]) > 1.0e-20)
	{
		v8 = a2[0];
		if (fabs(a2[1] * a2[1] + v8 * v8 + a2[2] * a2[2]) > 1.0e-20)
		{
			v9 = a1[1];
			v10 = a1[2];
			v11 = a2[2];
			v12 = a2[1];
			if (fabs(
				((v10 * v8 - v11 * v7) * (v10 * v8 - v11 * v7)
					+ (v11 * v9 - v12 * v10) * (v11 * v9 - v12 * v10)
					+ (v12 * v7 - v9 * v8) * (v12 * v7 - v9 * v8))
				/ ((v9 * v9 + v7 * v7 + v10 * v10)
					* (v12 * v12 + v8 * v8 + v11 * v11))) > a3 * a3)
				return 0;
		}
	}
	return 0;
}

static int UF_VEC3_ask_perpendicular1(const double a1[3], double a2[3])
{
	int result;
	double v4, v6, v7, v8;

	if (fabs(a1[0]) >= 1.0e19 || fabs(a1[1]) >= 1.0e19 || fabs(a1[2]) >= 1.0e19)
	{
		return 1;
	}
	v4 = a1[0];
	result = 0;
	if (fabs(a1[1] * a1[1] + v4 * v4 + a1[2] * a1[2]) <= 1.0e-20)
	{
		a2[0] = 0;
		a2[1] = 0;
		a2[2] = 0;
		return result;
	}
	v6 = fabs(v4);
	v7 = fabs(a1[1]);
	v8 = fabs(a1[2]) - 1.0e-10;
	if (v6 < v7 - 1.0e-10)
	{
		if (v6 < v8)
		{
			a2[0] = 0;
			if (v7 >= v8)
			{
				a2[1] = -a1[2];
				a2[2] = a1[1];
				return 1;
			}
			else
			{
				a2[1] = a1[2];
				result = 1;
				a2[2] = -a1[1];
			}
			return result;
		}
		a2[2] = 0;
		if (v7 >= v6 - 1.0e-10)
		{
			a2[0] = a1[1];
		LABEL_22:
			result = 1;
			a2[1] = -a1[0];
			return result;
		}
		goto LABEL_20;
	}
	if (v7 >= v8)
	{
		a2[2] = 0;
		if (v7 >= v6 - 1.0e-10)
		{
			a2[0] = a1[1];
			goto LABEL_22;
		}
	LABEL_20:
		a2[0] = -a1[1];
		a2[1] = a1[0];
		return 1;
	}
	a2[1] = 0;
	if (v6 >= v8)
	{
		a2[0] = -a1[2];
		a2[2] = a1[0];
		return 1;
	}
	else
	{
		a2[0] = a1[2];
		result = 1;
		a2[2] = -a1[0];
	}
	return result;
}

static int UF_VEC3_is_perpendicular1(const double a1[3], const double a2[3], double a3, int* result)
{
	int ret = 0;
	*result = 1;
	if (fabs(a1[0]) >= 1.0e19 || fabs(a1[1]) >= 1.0e19 || fabs(a1[2]) >= 1.0e19)
	{
		return 1;
	}
	if (fabs(a2[0]) >= 1.0e19 || fabs(a2[1]) >= 1.0e19 || fabs(a2[2]) >= 1.0e19)
	{
		return 1;
	}
	double t1 = a1[0] * a2[0] + a1[1] * a2[1] + a1[2] * a2[2];
	double t2 = a1[0] * a2[0] + a1[1] * a2[1] + a1[2] * a2[2];

	double t3 = a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2];
	double t4 = a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2];

	double v3 = a3 * a3;
	if (fabs(t1) > 1.0e-20 && fabs(t2) > 1.0e-20)
	{
		if (fabs(t1 * t2 / (t3 * t4)) > v3)
		{
			*result = 0;
		}
	}
	return ret;
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
	if (fabs(arc_coords.end_angle - arc_coords.start_angle - M_PI) < 0.001)
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

//create_flag=1= arc   2= circle
static int UF_CURVE_create_arc_thru_3pts(int     create_flag, double  first_point[3], double  second_point[3], double  third_point[3], int* arc_tag)
{
	szwPoint start = { first_point[0], first_point[1],  first_point[2] };
	szwPoint center = { second_point[0], second_point[1],  second_point[2] };
	szwPoint end = { third_point[0], third_point[1],  third_point[2] };
	int ret = 0;
	if (create_flag == 1)
	{
		ret = cvxPartArc3pt(&start,& end,&center, arc_tag);
	} 
	else if (create_flag == 2)
	{
		ret = cvxPartCir3pt(&start, &end, &center, arc_tag);
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

static int UF_MODL_ask_face_uv_minmax(int idFace, double(&uv_min_max)[4])
{
	svxLimit U, V;
	int ret = cvxFaceParam(idFace, &U, &V);
	uv_min_max[0] = U.min;
	uv_min_max[1] = U.max;
	uv_min_max[2] = V.min;
	uv_min_max[3] = V.max;
	return ret;
}

//radii没获取到
static int UF_MODL_ask_face_props(int face_id, double uv[2], double point[3], double u1[3], double v1[3], double u2[3], double v2[3], double unit_norm[3], double radii[2])
{
	radii[0] = 0;
	radii[1] = 0;
	/*
	svxPoint pt;
	svxVector dir;
	int iret = cvxFaceEval(face_id, uv[0], uv[1], &pt, &dir);
	*/
	szwEntityHandle faceHandle;
	id2hand(face_id, &faceHandle);
	szwSurface surfacedata;
	ZwFaceSurfaceDataGet(faceHandle, &surfacedata);

	szwPoint point1;
	szwVector normalDirection, uTangent, vTangent;
	int ret=ZwFacePointGetByUVparameter(surfacedata, uv[0], uv[1], &point1, &normalDirection,& uTangent, &vTangent);
	ZwSurfaceDataFree(&surfacedata);
	ZwEntityHandleFree(&faceHandle);

	point[0] = point1.x;
	point[1] = point1.y;
	point[2] = point1.z;

	unit_norm[0] = normalDirection.x;
	unit_norm[1] = normalDirection.y;
	unit_norm[2] = normalDirection.z;

	u1[0] = uTangent.x;
	u1[1] = uTangent.y;
	u1[2] = uTangent.z;

	v1[0] = vTangent.x;
	v1[1] = vTangent.y;
	v1[2] = vTangent.z;
	return ret;
}



/*
ZW_CURVE_LINE = 1
	ZW_CURVE_ARC = 2
	ZW_CURVE_CIRCLE = 3
	ZW_CURVE_NURB = 4
	ZW_CURVE_ELLIPSE2 = 5*/

//allcurve里面的还要循环释放	ZwEntityHandleFree(&entity);
static int UF_OBJ_cyl_curves(int type, vector <szwEntityHandle>& allcurve)
{
	allcurve.clear();
	int countCurve=0;
	szwEntityHandle* curveList=NULL;
	int ret=ZwCurveListGet(&countCurve, &curveList);
	for (size_t i = 0; i < countCurve; i++)
	{
		szwCurve curve;
		ZwCurveNURBSDataGet(curveList[i], 1, &curve);
		if (curve.type == type)
		{
			allcurve.push_back(curveList[i]);
		}
	}
	ZwEntityHandleListFree(countCurve, &curveList);
	return ret;
}

static int UF_OBJ_cycle_objs_in_part(int part_tag, int type, vector <int>& objs)
{
	int ret = 0;
	objs.clear();
	int count;
	int* objects;
	switch (type)
	{
	case UF_dummy_type:
	{
		ret = cvxPartInqAxis(&count, &objects);//平面，ZW没有
	}
	case UF_point_type:
	{
		ret = cvxPartInqPoints(&count, &objects);
	}
	case UF_line_type:
	{
		vector <szwEntityHandle> allcurve;
		UF_OBJ_cyl_curves(ZW_CURVE_LINE, allcurve);
		objs = hand2id(allcurve.size(), allcurve.data());
		for (size_t i = 0; i < allcurve.size(); i++)
		{
			ZwEntityHandleFree(&allcurve[i]);
		}
	}
	case UF_circle_type:
	{
		vector <szwEntityHandle> allcurve;
		UF_OBJ_cyl_curves(ZW_CURVE_ARC, allcurve);
		objs = hand2id(allcurve.size(), allcurve.data());
		for (size_t i = 0; i < allcurve.size(); i++)
		{
			ZwEntityHandleFree(&allcurve[i]);
		}
		UF_OBJ_cyl_curves(ZW_CURVE_CIRCLE, allcurve);
		vector <int> objs1;
		objs1 = hand2id(allcurve.size(), allcurve.data());
		for (size_t i = 0; i < allcurve.size(); i++)
		{
			objs.push_back(objs1[i]);
			ZwEntityHandleFree(&allcurve[i]);
		}
	}
	case UF_conic_type:
	{	
		vector <szwEntityHandle> allcurve;
		UF_OBJ_cyl_curves(ZW_CURVE_ELLIPSE2, allcurve);
		objs = hand2id(allcurve.size(), allcurve.data());
		for (size_t i = 0; i < allcurve.size(); i++)
		{
			ZwEntityHandleFree(&allcurve[i]);
		}
	}
	case UF_spline_type:
	{	
		/*
		ret = cvxPartInqCurves(&count, &objects);
		for (size_t i = 0; i < count; i++)
		{
			svxCurve Crv;
			cvxPartInqCurve(objects[i], 1, &Crv);
			if (VX_CRV_NURB == Crv.Type)
			{
				objs.push_back(objects[i]);
			}
		}*/
		vector <szwEntityHandle> allcurve;
		UF_OBJ_cyl_curves(ZW_CURVE_ELLIPSE2, allcurve);
		objs = hand2id(allcurve.size(), allcurve.data());
		for (size_t i = 0; i < allcurve.size(); i++)
		{
			ZwEntityHandleFree(&allcurve[i]);
		}
	}
	case UF_group_type:
	{
		ret = cvxPartInqGroupList(&count, &objects);
	}
	case UF_drafting_entity_type:
	{
		//cvxPartInqCurves(&count, &objects);
	}
	case UF_dimension_type:
	{
		//ZwDrawingSheetDimensionListGet
		//(int idDrawing, evxDimType * dimTypeList, int dimTypeCount, int* count, int** dims);
		//cvxDwgInqDims
		//cvxPartInqCurves(&count, &objects);
	}
	case UF_coordinate_system_type:
	{
		ret = cvxPartInqCsys(&count, &objects);
	}
	case UF_view_type:
	{
		ret = cvxPartInqViews(&count, &objects);
	}
	case UF_drawing_type:
	{
		int count;
		szwEntityHandle* sheetList;
		ZwDrawingSheetListGet(& count, &sheetList);
		objs = hand2id(count, sheetList);
		ZwEntityHandleListFree(count, &sheetList);
	}
	case UF_solid_type:
	{
		ret = cvxPartInqShapes(0, 0, &count, &objects);
	}

	case UF_sketch_type:
	{
		szwEntityHandle* sketchList = nullptr;
		ret = ZwSketchListGet(&count, &sketchList);
		for (size_t i = 0; i < count; i++)
		{
			int caotuid = hand2id(sketchList[i]);
		}
		ZwEntityHandleListFree(count, &sketchList);
	}
	case UF_texture_type:
	{
		ret = cvxPartInqTexts(&count, &objects);
	}
	case UF_feature_type:
	{
		szwEntityHandle* featureList = nullptr;
		ret = ZwFeatureListGet(&count, &featureList);
		for (size_t i = 0; i < count; i++)
		{
			int caotuid = hand2id(featureList[i]);
		}
		ZwEntityHandleListFree(count, &featureList);
	}

	default:
		break;
	}
	if (count>0)
	{
		for (size_t j = 0; j < count; j++)
		{
			objs.push_back(objects[j]);
		}
		cvxMemFree((void**)&objects);
	}
	return ret;
}

static int UF_OBJ_cycle_all(int part_tag,vector <int> & allobjects )
{

}

//这个参数与UG有一些差别
static int UF_CSYS_create_csys(double csys_origin[3],double matx[9], int * csysid)
{
	int idCSYS = 0;
	svxCSYSData CSYS = { 0 };
	cvxPartCSYSNewInit(VX_CSYS_ONLY_MATRIX, &CSYS);
	svxVector xdir = { matx[0], matx[1], matx[2] };
	svxVector ydir = { matx[3], matx[4], matx[5] };
	svxPoint pt = { csys_origin[0], csys_origin[1], csys_origin[2] };
	cvxMatPntVecs(&pt, &xdir, &ydir, &CSYS.Frame);
	int ret=cvxPartCSYSNew(&CSYS, csysid);
	return ret;
}

//这个参数与UG有一些差别
static int UF_CSYS_create_temp_csys(double csys_origin[3], double matx[9], int* csysid)
{
	int idCSYS = 0;
	svxCSYSData CSYS = { 0 };
	cvxPartCSYSNewInit(VX_CSYS_ONLY_MATRIX, &CSYS);
	svxVector xdir = { matx[0], matx[1], matx[2] };
	svxVector ydir = { matx[3], matx[4], matx[5] };
	svxPoint pt = { csys_origin[0], csys_origin[1], csys_origin[2] };
	cvxMatPntVecs(&pt, &xdir, &ydir, &CSYS.Frame);
	int ret = cvxPartCSYSNew(&CSYS, csysid);
	return ret;
}


static int UF_CSYS_set_wcs(int csysid)
{
	int ret = cvxPartActiveAsLCS(csysid);
	return ret;
}

//这个参数与UG有一些差别
static void UF_CSYS_ask_wcs(svxMatrix Mat)
{
	 cvxPartInqLCSMat(&Mat);
}

//这个参数与UG有一些差别
static int UF_CSYS_create_matrix(double csys_origin[3],double matx[9], svxMatrix& Mat)
{
	svxVector xdir = { matx[0], matx[1], matx[2] };
	svxVector ydir = { matx[3], matx[4], matx[5] };
	svxPoint pt = { csys_origin[0], csys_origin[1], csys_origin[2] };
	int ret = cvxMatPntVecs(&pt, &xdir, &ydir, &Mat);
	return ret;
}

//这个参数与UG有一些差别
static int UF_CSYS_create_matrix(double matx[9], svxMatrix& Mat)
{
	svxVector xdir = { matx[0], matx[1], matx[2] };
	svxVector ydir = { matx[3], matx[4], matx[5] };
	svxPoint pt = { 0};
	int ret = cvxMatPntVecs(&pt, &xdir, &ydir, &Mat);
	return ret;
}

//这个参数与UG有一些差别
static int UF_CSYS_ask_matrix_values(svxMatrix matx,double matx_vec[9])
{
	matx_vec[0] = matx.xx;
	matx_vec[1] = matx.xy;
	matx_vec[2] = matx.xz;

	matx_vec[3] = matx.yx;
	matx_vec[4] = matx.yy;
	matx_vec[5] = matx.yz;

	matx_vec[6] = matx.zx;
	matx_vec[7] = matx.zy;
	matx_vec[8] = matx.zz;
	return 0;
}

//这个只是查面的
static int  UF_CSYS_ask_matrix_of_object(int obj,double matx_vec[9])
{
	svxSrfPrim SrfPrim;
	cvxPartInqFaceSrfPrim(obj, &SrfPrim);
	if (SrfPrim.srfType == VX_SF_PRIM_PLN)
	{
		UF_CSYS_ask_matrix_values(SrfPrim.srfData.pln.form, matx_vec);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_SPH)
	{
		UF_CSYS_ask_matrix_values(SrfPrim.srfData.sph.form, matx_vec);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_CON)
	{
		UF_CSYS_ask_matrix_values(SrfPrim.srfData.con.form, matx_vec);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_CYL)
	{
		UF_CSYS_ask_matrix_values(SrfPrim.srfData.cyl.form, matx_vec);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_NURB)
	{
		svxMatrix Mat;
		cvxMatInit(&Mat);
		UF_CSYS_ask_matrix_values(Mat, matx_vec);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_ELLSO)
	{
		UF_CSYS_ask_matrix_values(SrfPrim.srfData.ellso.form, matx_vec);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_TORUS)
	{
		UF_CSYS_ask_matrix_values(SrfPrim.srfData.torus.form, matx_vec);
	}
	cvxSrfPrimFree(&SrfPrim);
	return 0;
}



static int UF_CSYS_ask_csys_info(int csys_id, double matrix[9], double csys_origin[3])
{
	svxCSYSData CSYS;
	int ret =cvxCSYSGetData(csys_id, &CSYS);
	UF_CSYS_ask_matrix_values(CSYS.Frame, matrix);
	csys_origin[0] = CSYS.Frame.xt;
	csys_origin[1] = CSYS.Frame.yt;
	csys_origin[2] = CSYS.Frame.zt;
	return ret;
}



//圆锥点是中间的，直径也是中间的，ZW提供是最大的
static	int UF_MODL_ask_face_data(int face_id, int* facetype, double point[], double dir[], double box[], double* radius, double* rad_data, int* norm_dir)
{
	*facetype = 0;
	for (size_t i = 0; i < 3; i++)
	{
		point[i] = 0;
		dir[i] = 0;
		box[i] = 0;
		box[i + 3] = 0;
	}
	//ZW-0是凸面，1是凹面
	int aotu = cvxFaceIsConcave(face_id);
	if (aotu == 0)
	{
		*norm_dir = 1;
	}
	else if (aotu == 1)
	{
		*norm_dir = 0;
	}
	svxBndBox Box;
	svxSrfPrim SrfPrim;
	cvxPartInqFaceSrfPrim(face_id, &SrfPrim);
	if (SrfPrim.srfType == VX_SF_PRIM_PLN)
	{
		*facetype = UF_MODL_PLANAR_FACE;
		UF_MTX3_point(SrfPrim.srfData.pln.form, point);
		UF_MTX3_z_vec(SrfPrim.srfData.pln.form, dir);
		cvxPartInqEntBox(face_id, &SrfPrim.srfData.pln.form, &Box);
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_SPH)
	{
		*facetype = UF_MODL_SPHERICAL_FACE;
		UF_MTX3_point(SrfPrim.srfData.sph.form, point);
		UF_MTX3_z_vec(SrfPrim.srfData.sph.form, dir);
		*radius = SrfPrim.srfData.sph.radius;
		//cvxPartInqEntBox(face_id, &SrfPrim.srfData.sph.form, &Box);//与NX一样吧，不按圆柱方向获取大小了，按绝对坐标取
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_CON)
	{
		*facetype = UF_MODL_CONICAL_FACE;
		svxBndBox conBox;
		cvxPartInqEntBox(face_id, &SrfPrim.srfData.con.form, &conBox);
		double len = conBox.Z.max - conBox.Z.min;
		svxPoint cenpt = { SrfPrim.srfData.con.form.xt,SrfPrim.srfData.con.form.yt,SrfPrim.srfData.con.form.zt };
		svxVector zvec = { SrfPrim.srfData.con.form.zx,SrfPrim.srfData.con.form.zy,SrfPrim.srfData.con.form.zz };
		UF_MTX3_z_vec(SrfPrim.srfData.con.form, dir);
		cvxPntTranslate(&cenpt, &zvec, len/2);
		UF_VEC3_copy(cenpt, point);
		double duanlen = SrfPrim.srfData.con.radius2-SrfPrim.srfData.con.radius1;
		double xielen = sqrt(duanlen * duanlen + len * len);
		double angle = std::atan(duanlen/len);
		double degrees = angle * (180.0 / M_PI);
		double zlen = fabs(tan(angle) * len / 2);

		*radius = SrfPrim.srfData.con.radius2 - zlen;//大半径跟UG不一样，UG要
		*rad_data = angle;//小半径//要转角度

		//char msg[256] = "";
		//sprintf_s(msg, "r1=%f r2=%f r0=%f  zlen=%f ", SrfPrim.srfData.con.radius1, SrfPrim.srfData.con.radius2, SrfPrim.srfData.con.radius2- zlen, zlen);
		//cvxMsgDisp(msg);

		//sprintf_s(msg, "r1=%f r2=%f dl=%f l=%f xl=%f a=%f j=%f", SrfPrim.srfData.con.radius1, SrfPrim.srfData.con.radius2, duanlen, len, xielen, angle, degrees);
		//cvxMsgDisp(msg);
		//长度有问题
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_CYL)
	{
		*facetype = UF_MODL_CYLINDRICAL_FACE;
		svxBndBox cyl;
		cvxPartInqEntBox(face_id, &SrfPrim.srfData.cyl.form, &cyl);
		double len = cyl.Z.max - cyl.Z.min;
		svxPoint cenpt = { SrfPrim.srfData.con.form.xt,SrfPrim.srfData.cyl.form.yt,SrfPrim.srfData.con.form.zt };
		svxVector zvec = { SrfPrim.srfData.con.form.zx,SrfPrim.srfData.cyl.form.zy,SrfPrim.srfData.con.form.zz };
		cvxPntTranslate(&cenpt, &zvec, len / 2);
		UF_VEC3_copy(cenpt, point);
		UF_MTX3_z_vec(SrfPrim.srfData.cyl.form, dir);
		*radius = SrfPrim.srfData.cyl.radius;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_NURB)
	{
		*facetype = UF_MODL_PARAMETRIC_FACE;
		//SrfPrim.srfData.srf.P.
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_ELLSO)
	{
		*facetype = UF_MODL_FOREIGN_FACE;
		UF_MTX3_point(SrfPrim.srfData.ellso.form, point);
		UF_MTX3_z_vec(SrfPrim.srfData.ellso.form, dir);
		*radius = SrfPrim.srfData.ellso.xlen;
		*rad_data = SrfPrim.srfData.ellso.ylen;
	}
	else if (SrfPrim.srfType == VX_SF_PRIM_TORUS)
	{
		//圆环面
		*facetype = UF_MODL_TOROIDAL_FACE;
		UF_MTX3_point(SrfPrim.srfData.torus.form, point);
		UF_MTX3_z_vec(SrfPrim.srfData.torus.form, dir);
		*radius = SrfPrim.srfData.torus.dPathRadius;
		*rad_data = SrfPrim.srfData.torus.dProfRadius;
	}
	cvxSrfPrimFree(&SrfPrim);

	svxMatrix Mat;
	cvxMatInit(&Mat);
	cvxPartInqEntBox(face_id, &Mat, &Box);
	box[0] = Box.X.min;
	box[1] = Box.Y.min;
	box[2] = Box.Z.min;
	box[3] = Box.X.max;
	box[4] = Box.Y.max;
	box[5] = Box.Z.max;
	return 0;
}

//体找特征
static int UF_MODL_ask_body_features(int idshape,vector <int> &features )
{
	int Count = 0;
	int* Features = NULL;
	int ret =cvxPartInqShapeFtrs(idshape, 0, &Count,& Features);
	for (size_t i = 0; i < Count; i++)
	{
		features.push_back(Features[i]);
	}
	if (Features!=NULL)
	{
		cvxMemFree((void**)&Features);
	}
	return ret;
}

//特征找体
static int UF_MODL_ask_feat_body(int Features,int *idshape)
{
	int* label = NULL;
	cvxEntLabelGet(Features, &label);
	cvxEntByLabel(label, 1, idshape);
	cvxMemFree((void**)&label);
	return 0;
}

//特征找面
static int UF_MODL_ask_feat_faces(int feat_id, vector <int>& faces)
{
	int n_ents = 0;
	int* ents_ids = NULL;
	int ret = cvxPartInqFtrEnts(feat_id, VX_ENT_FACE, &n_ents, &ents_ids);
	for (size_t j = 0; j < n_ents; j++)
	{
		faces.push_back(ents_ids[j]);
	}
	if (ents_ids != NULL)
	{
		cvxMemFree((void**)&ents_ids);
	}
	return ret;
}

//特征找边
static int UF_MODL_ask_feat_edges(int feat_id, vector <int>& edges)
{
	int n_ents = 0;
	int* ents_ids = NULL;
	int ret = cvxPartInqFtrEnts(feat_id, VX_ENT_EDGE, &n_ents, &ents_ids);
	for (size_t j = 0; j < n_ents; j++)
	{
		edges.push_back(ents_ids[j]);
	}
	if (ents_ids != NULL)
	{
		cvxMemFree((void**)&ents_ids);
	}
	return ret;
}


//特征找所有体
static int UF_MODL_ask_feat_bodys(int feat_id, vector <int>& bodys)
{
	int n_ents = 0;
	int* ents_ids = NULL;
	int ret = cvxPartInqFtrEnts(feat_id, VX_ENT_SHAPE, &n_ents, &ents_ids);
	for (size_t j = 0; j < n_ents; j++)
	{
		bodys.push_back(ents_ids[j]);
	}
	if (ents_ids != NULL)
	{
		cvxMemFree((void**)&ents_ids);
	}
	return ret;
}

//多选
static int UF_UI_select_with_class_dialog(char* message, char* title, int scope, evxEntInpOpt option, void* user_data, int* response, int* count, int** object)
{
	int ret = cvxGetEnts(title, option, 1, count, object);
	return ret;
}

//单选择
static int UF_UI_select_with_single_dialog(char* message, char* title, int scope, evxEntInpOpt option, void* user_data, int* response, int* object, double cursor[3], int* view)
{
	cursor[0] = 0;
	cursor[1] = 0;
	cursor[2] = 0;
	*view = 0;
	int ret = cvxGetEnt(title, option, 1, object);
	return ret;
}


//属性删除
static int UF_ATTR_delete(int object, int type, char* title)
{

}

//与UG有一点点差别
//获取实体属性
// 注意：如果是字符串使用完毕后应该释放内存	// free(values[i].value.string);
static int UF_ATTR_cycle(int obj_tag, int* indx, int type, char* title,vector<UF_ATTR_value_t> &values)
{
	indx = 0;	
	//获取所有属性
	svxPartAttribute At;
	cvxShellAtGet(obj_tag, &At);
	for (int k = 0; k < At.UserAttributeCount; k++)
	{
		if (strcmp(At.UserAttribute[k].label, title)==0)
		{
			UF_ATTR_value_t att;

			if (At.UserAttribute[k].type == VX_ATTR_INT)
			{
				att.type = UF_ATTR_integer;
				att.value.integer = (int)At.UserAttribute[k].dValue;
			}
			else if (At.UserAttribute[k].type == VX_ATTR_REAL)
			{
				att.type = UF_ATTR_real;
				att.value.real = At.UserAttribute[k].dValue;
			}
			else if (At.UserAttribute[k].type == VX_ATTR_STRING)
			{
				att.type = UF_ATTR_string;
				att.value.string = (char*)malloc(UF_ATTR_MAX_STRING_LEN * sizeof(char));
				strcpy_s(att.value.string, UF_ATTR_MAX_STRING_LEN, At.UserAttribute[k].strValue);
				// 注意：使用完毕后应该释放内存	// free(str);
			}
			else if (At.UserAttribute[k].type == VX_ATTR_DATE)
			{
				att.type = UF_ATTR_time;
				att.value.string = (char*)malloc(UF_ATTR_MAX_STRING_LEN * sizeof(char));
				strcpy_s(att.value.string, UF_ATTR_MAX_STRING_LEN,At.UserAttribute[k].strValue);
				// 注意：使用完毕后应该释放内存	// free(str);
			}
			
			values.push_back(att);
		}
	}

	return 0;
}

//设置字符串属性
static void set_obj_attr(int object, const char* title, const char* vlaue)
{
	evxAtItemId itemId = VX_AT_USER;
	svxAttribute At;
	At.type = VX_ATTR_STRING;
	strcpy_s(At.label, title);
	strcpy_s(At.strValue, vlaue);
	cvxShellAtItemSet(object, itemId, &At);
}

//设置int属性
static void set_obj_attr(int object, const char* title, int vlaue)
{
	evxAtItemId itemId = VX_AT_USER;
	svxAttribute At;
	At.type = VX_ATTR_STRING;
	strcpy_s(At.label, title);
	At.dValue= vlaue;
	cvxShellAtItemSet(object, itemId, &At);
}

//设置double属性
static void set_obj_attr(int object, const char* title, double vlaue)
{
	evxAtItemId itemId = VX_AT_USER;
	svxAttribute At;
	At.type = VX_ATTR_STRING;
	strcpy_s(At.label, title);
	At.dValue = vlaue;
	cvxShellAtItemSet(object, itemId, &At);
}

//设置属性
static int UF_ATTR_assign(int object, const char* title, UF_ATTR_value_t vlaue)
{
	if (vlaue.type== UF_ATTR_integer)
	{
		set_obj_attr(object, title, vlaue.value.integer);
	}
	else if (vlaue.type == UF_ATTR_real)
	{
		set_obj_attr(object, title, vlaue.value.real);
	}
	else if (vlaue.type == UF_ATTR_string)
	{
		set_obj_attr(object, title, vlaue.value.string);
	}
	return 0;
}


struct UF_MODL_ray_hit_point_info_s
{
	double    hit_point[3];
	int       hit_face;
};

typedef struct UF_MODL_ray_hit_point_info_s UF_MODL_ray_hit_point_info_t,
* UF_MODL_ray_hit_point_info_p_t;

//射线
static int UF_MODL_trace_a_ray(int num_bodies, int* bodys, double pt[3], double dir[3], double trans[16], int num_desired, int* num_results, vector <UF_MODL_ray_hit_point_info_t>& hit_list)
{
	szwEntityHandle* shapes;
	ZwEntityIdTransfer(num_bodies, bodys, shapes);
	ezwFaceTrim faceTrim = ZW_TRIM_ALL;
	szwAxis ray;
	ray.point = { pt[0],  pt[1],  pt[2] };
	int pttag;
	ray.direction = { dir[0], dir[1], dir[2] };
	int infinite = 0;
	double length = 0.0;
	int count;
	szwIntersectionPoint* intersectionPoints;
	int ret = ZwRayShapeIntersect(num_bodies, shapes, faceTrim, ray, infinite, length, &count, &intersectionPoints);

	if (count > 0)
	{
		*num_results = count;
		for (size_t i = 0; i < count; i++)
		{
			int face = hand2id(intersectionPoints[i].faceHandle);
			UF_MODL_ray_hit_point_info_t hitobj;
			hitobj.hit_face = face;
			hitobj.hit_point[0] = intersectionPoints[i].point.x;
			hitobj.hit_point[1] = intersectionPoints[i].point.y;
			hitobj.hit_point[2] = intersectionPoints[i].point.z;
			hit_list.push_back(hitobj);
		}
		ZwMemoryFree((void**)&intersectionPoints);
		ZwEntityHandleListFree(num_bodies, &shapes);
	}
	return  ret;
}

//射线
static int UF_MODL_trace_a_ray(int num_bodies, szwEntityHandle* shapes, double pt[3], double dir[3], vector <UF_MODL_ray_hit_point_info_t>& hit_list)
{
	ezwFaceTrim faceTrim = ZW_TRIM_ALL;
	szwAxis ray;
	ray.point = { pt[0],  pt[1],  pt[2] };
	int pttag;
	ray.direction = { dir[0], dir[1], dir[2] };
	int infinite = 0;
	double length = 0.0;
	int count;
	szwIntersectionPoint* intersectionPoints;
	int ret = ZwRayShapeIntersect(num_bodies, shapes, faceTrim, ray, infinite, length, &count, &intersectionPoints);

	if (count > 0)
	{
		for (size_t i = 0; i < count; i++)
		{
			int face = hand2id(intersectionPoints[i].faceHandle);
			UF_MODL_ray_hit_point_info_t hitobj;
			hitobj.hit_face = face;
			hitobj.hit_point[0] = intersectionPoints[i].point.x;
			hitobj.hit_point[1] = intersectionPoints[i].point.y;
			hitobj.hit_point[2] = intersectionPoints[i].point.z;
			hit_list.push_back(hitobj);
		}
		ZwMemoryFree((void**)&intersectionPoints);
		ZwEntityHandleListFree(num_bodies, &shapes);
	}
	return  ret;
}

//射线
static int UF_MODL_trace_a_ray(int num_bodies, int* bodys, double pt[3], double dir[3], vector <UF_MODL_ray_hit_point_info_t>& hit_list)
{
	szwEntityHandle* shapes;
	ZwEntityIdTransfer(num_bodies, bodys, shapes);
	int ret = UF_MODL_trace_a_ray(num_bodies, shapes, pt, dir, hit_list);
	return  ret;
}

//按方向矩阵查包围框
static int UF_MODL_ask_bounding_box_exact(int idShape, int csys_id, double  min_corner[3], double  directions[3][3], double  distances[3])
{
	svxCSYSData csys_data;
	cvxCSYSGetData(csys_id, &csys_data);
	svxMatrix matrix = csys_data.Frame;
	svxVector vec[6];
	UF_MTX3_x_vec(matrix, vec[0]);
	UF_MTX3_y_vec(matrix, vec[1]);
	UF_MTX3_z_vec(matrix, vec[2]);
	for (size_t i = 3; i < 6; i++)
	{
		vec[i] = vec[i - 3];
		cvxVecReverse(&vec[i]);
	}
	svxMatrix mat;
	cvxMatInit(&mat);
	svxPoint boxpt[6];
	svxPoint boxwcspt[6];
	for (size_t i = 0; i < 6; i++)
	{
		int Count;
		int* idEnts;
		svxPoint Point;
		cvxPartInqShapeExtreme(idShape, &vec[i], &Count, &idEnts, &boxpt[i]);
		cvxMemFree((void**)&idEnts);
		UF_CSYS_map_point(mat, boxpt[i], matrix, boxwcspt[i]);
		//int pttag;
		//cvxPartPnt(&Point, &pttag);
	}
	double minX = boxwcspt[3].x;
	double minY = boxwcspt[4].y;
	double minZ = boxwcspt[5].z;
	double maxX = boxwcspt[0].x;
	double maxY = boxwcspt[1].y;
	double maxZ = boxwcspt[2].z;

	distances[0] = maxX - minX;
	distances[1] = maxY - minY;
	distances[2] = maxZ - minZ;

	svxPoint minipt = { minX ,minY,minZ }, abspt;
	UF_CSYS_map_point(matrix, minipt, mat, abspt);

	min_corner[0] = abspt.x;
	min_corner[1] = abspt.y;
	min_corner[2] = abspt.z;

	return 0;
}

//按方向矩阵查包围框
static int UF_MODL_ask_bounding_box_exact(int idShape, svxMatrix matrix, double  min_corner[3], double  cenpt[3], double  distances[3])
{
	svxBndBox box;
	cvxPartInqEntBox(idShape, &matrix, &box);
	svxVector vec[6];
	UF_MTX3_x_vec(matrix, vec[0]);
	UF_MTX3_y_vec(matrix, vec[1]);
	UF_MTX3_z_vec(matrix, vec[2]);
	for (size_t i = 3; i < 6; i++)
	{
		vec[i] = vec[i - 3];
		cvxVecReverse(&vec[i]);
	}
	svxMatrix mat;
	cvxMatInit(&mat);
	svxPoint boxpt[6];
	svxPoint boxwcspt[6];
	for (size_t i = 0; i < 6; i++)
	{
		int Count;
		int* idEnts;
		svxPoint Point;
		cvxPartInqShapeExtreme(idShape, &vec[i], &Count, &idEnts, &boxpt[i]);
		cvxMemFree((void**)&idEnts);
		UF_CSYS_map_point(mat, boxpt[i], matrix, boxwcspt[i]);
		//int pttag;
		//cvxPartPnt(&Point, &pttag);
	}
	double minX = boxwcspt[3].x;
	double minY = boxwcspt[4].y;
	double minZ = boxwcspt[5].z;
	double maxX = boxwcspt[0].x;
	double maxY = boxwcspt[1].y;
	double maxZ = boxwcspt[2].z;

	distances[0] = maxX - minX;
	distances[1] = maxY - minY;
	distances[2] = maxZ - minZ;

	svxPoint minipt = { minX ,minY,minZ }, abspt;
	UF_CSYS_map_point(matrix, minipt, mat, abspt);

	min_corner[0] = abspt.x;
	min_corner[1] = abspt.y;
	min_corner[2] = abspt.z;

	svxPoint wcscenpt = { minX + distances[0] / 2 ,minY + distances[1] / 2,minZ + distances[2] / 2 };
	UF_CSYS_map_point(matrix, wcscenpt, mat, abspt);

	cenpt[0] = abspt.x;
	cenpt[1] = abspt.y;
	cenpt[2] = abspt.z;

	return 0;
}


//查包围框
static int UF_MODL_ask_bounding_box(int idShape, double bounding_box[6])
{
	svxBndBox Box1 = { 0 };
	int ret = cvxEntBndBox(idShape, &Box1);
	bounding_box[0] = Box1.X.min;
	bounding_box[1] = Box1.Y.min;
	bounding_box[2] = Box1.Z.min;

	bounding_box[3] = Box1.X.max;
	bounding_box[4] = Box1.Y.max;
	bounding_box[5] = Box1.Z.max;
	return ret;
}

//边倒圆角
static int UF_MODL_create_edge_blend(int obj, double rad)
{
	int ret=cvxPartFillet(1, &obj, rad);
	return ret;
}

//多条边倒圆角
static int UF_MODL_create_edge_blend(int count ,int *objs, double rad)
{
	int ret = cvxPartFillet(count, objs, rad);
	return ret;
}
//边倒圆角
static int  UF_MODL_create_blend(int obj, double rad)
{
	int ret = cvxPartFillet(1, &obj, rad);
	return ret;
}
//边倒圆角
static int  UF_MODL_create_blend1(int obj, double rad)
{
	int ret = cvxPartFillet(1, &obj, rad);
	return ret;
}

//倒斜角
static int UF_MODL_create_chamfer(int subtype, char* offset1, char* offset2, char* theta, vector <int> edges, tag_t* feature_obj_id)
{
	int ret = 1;
	double ang = atof(theta);
	if (fabs(ang-45)<0.01 && subtype==1)
	{
		double  chamfer_size = atof(offset1);
		ret=cvxPartChamConst(edges.size(), edges.data(), chamfer_size);
	}
	return ret;
}

//查询透明度
static int UF_OBJ_ask_translucency(int face, int *value)
{
	svxFaceAt At;
	int ret = cvxPartInqFaceAt(face, &At);
	*value = At.trans;
	return ret;
}

//设置透明度
static int UF_OBJ_set_translucency(int face, int value)
{
	int ret = cvxEntTransSet(value, 1, &face);
	return ret;
}

static int UF_MODL_unite_bodies(int target, int tool)
{
	int ret = cvxPartBool(VX_BOOL_ADD, target, 1, &tool, 0);
	return ret;
}

static int UF_MODL_subtract_bodies(int  target, int  tool, int* num_result, int** resulting_bodies)
{
	int opstart = cvxOpCount();
	int ret = cvxPartBool(VX_BOOL_REMOVE, target, 1, &tool, 0);
	cvxEntNewAll(opstart, VX_ENT_SHAPE, num_result, resulting_bodies);
	if (*num_result == 0)
	{
		*resulting_bodies[0] = target;
	}
	return ret;
}

static int UF_MODL_intersect_bodies(int  target, int  tool, int* num_result, int** resulting_bodies)
{
	int opstart = cvxOpCount();
	int ret = cvxPartBool(VX_BOOL_INTERSECT, target, 1, &tool, 0);
	cvxEntNewAll(opstart, VX_ENT_SHAPE, num_result, resulting_bodies);
	if (*num_result==0)
	{
		*resulting_bodies[0] = target;
	}
	return ret;
}

static int UF_MODL_unite_bodies_with_retained_options(int target, int tool, logical retain_target_body, logical  retain_tool_body, int* frec_eid)
{
	int keep = 0;
	if (retain_tool_body)
	{
		keep = 1;
	}
	int Count, Count1;
	int* Feats1, * Feats;
	cvxPartInqShapeFtrs(target, 0, &Count, &Feats);
	int ret = cvxPartBool(VX_BOOL_ADD, target, 1, &tool, keep);
	cvxPartInqShapeFtrs(target, 0, &Count1, &Feats1);	if (Count1 > Count)
	{
		for (size_t i = 0; i < Count1; i++)
		{
			bool isnewfeat = true;
			for (size_t n = 0; n < Count; n++)
			{
				if (Feats1[i] == Feats[n])
				{
					isnewfeat = false;
					break;
				}
			}
			if (isnewfeat)
			{
				*frec_eid = Feats1[i];
				break;
			}
		}
	}
	return ret;
}

static int UF_MODL_subtract_bodies_with_retained_options(int target,int tool, logical retain_target_body,logical  retain_tool_body, int* frec_eid)
{
	int keep = 0;
	if (retain_tool_body)
	{
		keep = 1;
	}
	int Count, Count1;
	int* Feats1, * Feats;
	cvxPartInqShapeFtrs(target, 0, &Count, &Feats);
	int ret = cvxPartBool(VX_BOOL_REMOVE, target, 1, &tool, keep);
	cvxPartInqShapeFtrs(target, 0, &Count1, &Feats1);
	if (Count1 > Count)
	{
		for (size_t i = 0; i < Count1; i++)
		{
			bool isnewfeat = true;
			for (size_t n = 0; n < Count; n++)
			{
				if (Feats1[i] == Feats[n])
				{
					isnewfeat = false;
					break;
				}
			}
			if (isnewfeat)
			{
				*frec_eid = Feats1[i];
				break;
			}
		}
	}
	return ret;
}

static int UF_MODL_intersect_bodies_with_retained_options(int target, int tool, logical retain_target_body, logical  retain_tool_body, int* frec_eid)
{
	int keep = 0;
	if (retain_tool_body)
	{
		keep = 1;
	}
	int Count, Count1;
	int* Feats1,*Feats;
	cvxPartInqShapeFtrs(target, 0, &Count,&Feats);
	int ret = cvxPartBool(VX_BOOL_INTERSECT, target, 1, &tool, keep);
	cvxPartInqShapeFtrs(target, 0, &Count1, &Feats1);
	if (Count1 > Count)
	{
		for (size_t i = 0; i < Count1; i++)
		{
			bool isnewfeat = true;
			for (size_t n = 0; n < Count; n++)
			{
				if (Feats1[i] == Feats[n])
				{
					isnewfeat = false;
					break;
				}
			}
			if (isnewfeat)
			{
				*frec_eid = Feats1[i];
				break;
			}
		}
	}
	return ret;
}


static int UF_MODL_ask_point_containment(double pt[3],int target, int* pt_status)
{
	svxPoint Pnt = { pt[0],pt[1], pt[2] };
	*pt_status = cvxPntIsOn(&Pnt, target);
	return *pt_status;
}