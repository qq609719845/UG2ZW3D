
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

//PK_BODY_unite_bodies(pk_cyl1, 1, booer_body, &n_bodies, &bodies);
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
