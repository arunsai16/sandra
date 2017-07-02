
/********************************************************************** ********************************
*                                         Shell program for CEG 4230/6230 EE/ME 4560/6560              * 
*                                         TO BE USED FROM LAB2 ONWARDS!                                * 
*   This is the basic shell you will need to program for the lab. Feel free to modify it in any        * 
*  way you desire, but remember it must maintain the menu driven format. You will need  to             * 
*  determine what each subroutine requires as far as parameter passing is concerned. For               * 
*  the moment, do not modify the NEST() or RESET() routines. These will be modified                    * 
*  later if required. Feel free to ask me any questions on this shell.                                 * 
********************************************************************************************************/ 

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

#include	"RobotSerialInterface.h"

int RobotMenu( void );
 
void GET_ANGLES ( void ); 
void MOD_ANGLES (void); 

void MOVE_ROBOT ( void ); 
void CALC_TRANSFORMS (double, double, double, double, double); 
void PRINT_POS(void); 
void printTransforms( void ); 
void INVERSE_KINEMATICS ( void ); 
void GRASP_OBJECT( void ); 

void CRWait( void );
double toRadians(double ang); // changed to double
double toDegrees(double rad); // used in printTransforms()

// cosine and sine shorthand functions
double C(double degrees);
double S(double degrees);

// Global variables

int numPostions = 0;     // number of positions
float commands[100][8];  // commands array
float curr[5] = {0,0,0,0,0};
double t5[5][5];
double t6[5][5];

int main( void ) 
{ 
    //double Pi = 4*arctan(1);
	int x;
    int option;
    int done;
    //int Status;
    unsigned int    RobotResult;
    //float   TempFloat;
    unsigned int    OutputOptions;
 
    printf("entering main program\n");
    done = 0; 
    
    OutputOptions = RVM1_WriteToFile | RVM1_WriteToComm | RVM1_WriteToScreen;
    RobotResult = RVM1_InitializeComm( OutputOptions , 
                                        "c:\\_rvm1_logs\\WEDTest.txt" );

    RVM1_FlushBuffers();
/**    RVM1_Reset();**/
	RVM1_Speed( 4 );
	RVM1_GripClose();
    RVM1_Nest();

    while (!done) 
    { 
        option = RobotMenu(); 
        switch (option) 
        { 
            case 1: 
                RVM1_Nest();
                break; 

            case 2: 
                GET_ANGLES(); 
                break; 

            case 3: 
                MOD_ANGLES(); 
                break; 

            case 4: 
                PRINT_POS(); 
		break; 

            case 5: 
		printTransforms(); 
                break; 

            case 6: 
		MOVE_ROBOT ();
                break; 

            case 7: 
                RVM1_Reset();
                RVM1_FlushBuffers();
                RVM1_Reset();
                RVM1_Origin();
				for(x=0; x<5; x++)
				{
					curr[x]=0;
				}
               
                break; 
            
            case 8:
				{ 
                INVERSE_KINEMATICS(); 
                break;
				} 
            
            case 9: 
                {
		GRASP_OBJECT(); 
                break; 
		} 
            case 10:
                 {
                       RVM1_Origin();
                       break;
                       }

            case 11: 
                 //RVM1_Nest();
                done = 1;
                break; 
                
            default:
                break; 
        }
    }

    RVM1_Nest();
    RVM1_FlushBuffers();
/**    RVM1_Reset();**/

    printf ("Exiting program.\n"); 
    exit(0); 
} 


/**********************************************************************/ 
int RobotMenu( void ) 
/* PRINTS USERS OPTIONS AND RETURNS THAT OPTION TO MAIN() */
{ 
    int OptionSelect;
    int DoneSelecting;
    DoneSelecting = 0; 
    while( !DoneSelecting ) 
    { 
        system("cls"); 
        printf ("\n\n"); 
        printf ("\t\t\t\tMAIN MENU\n"); 
        printf ("\t\t\t\t OPTIONS\n"); 
        printf ("\t\t1) Nest Robot\n"); 
        printf ("\t\t2) Input Positions\n"); 
        printf ("\t\t3) Modify Positions\n"); 
        printf ("\t\t4) Display Positions\n"); 
        printf ("\t\t5) Print Transforms\n"); 
        printf ("\t\t6) Move Robot to Stored Position(s)\n"); 
        printf ("\t\t7) Reset Robot\n"); 
        printf ("\t\t8) Inverse Kinematics\n");
		printf ("\t\t9) Grasp Object\n"); 
		printf ("\t\t10) Go To Origin\n"); 
        printf ("\t\t11) Quit Program.\n"); 
        printf ("\n"); 
        printf ("\t\t\tPlease enter option: "); 
        scanf("%d", &OptionSelect); 

        if (OptionSelect < 1 || OptionSelect > 11) 
        { 
            printf ("\n\n\t\t\tError number must be between 1"); 
            printf (" and %d",11); 
        } 
        else
        {
            DoneSelecting = 1;
        }
        printf ("\n\n\n"); 
    } 
    return( OptionSelect ); 
}


/**********************************************************************/ 
void CRWait(void)
{
	int		InputChar;

	fflush(stdin);
	printf("\t\tHit RETURN to Continue -> ");
	InputChar = getchar();
	fflush(stdin);
}

/**********************************************************************/ 

void INVERSE_KINEMATICS ()
{ 
	
	getchar();
}

/**************************************************************/
void GET_ANGLES() 
{
/* This subroutine prompts the user to enter in position data. Your code 
            should be written so that at least five positions will be stored. 
            As each angle is entered, verify it is within the joint range. In 
            addition, you should call CALC_TRANSFORMS() after each position is 
            entered and verify the robot will not crash into the table. Remember 
            you need only check the z position on both T05 and T06 to ensure crash protection. */ 
	float t;   //temp
	int good;
	int x, i=0;
	double x2, x3, x4;  // changed to double
	double z5, z6;      // changed to double

	printf("Enter number of positions you wish to enter: ");
	scanf("%d", &i);

	for(x=numPostions; x<(i+numPostions); x++)
	{
		good=0;
		while(!good)
		{
			printf("Enter base angle (-150 to 150): ");
			scanf("%f", &t);
			if (-150<=t && t<=150)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][0]=t;		
		
		//ToDo:Shoulder Angle
		good=0;
		while(!good)
		{
			printf("Enter shoulder angle (-30 to 105): ");
			scanf("%f", &t);
			if (-30<=t && t<=105)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][1]=t;	
		//ToDo:Elbow Angle
		good=0;
		while(!good)
		{
			printf("Enter Elbow angle (-105 to 2): ");
			scanf("%f", &t);
			if (-105<=t && t<=2)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][2]=t;	
		//ToDo:Wrist Pitch
		good=0;
		while(!good)
		{
			printf("Enter Wrist Pitch angle (-90 to 90): ");
			scanf("%f", &t);
			if (-90<=t && t<=90)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][3]=t;	
		//ToDo:Wrist Roll
		good=0;
		while(!good)
		{
			printf("Enter Wrist Roll angle (-180 to 180): ");
			scanf("%f", &t);
			if (-180<=t && t<=180)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][4]=t;	
		
		good=0;
		while(!good)
		{
			printf("Gripper before move(Open = 1, Close = -1): ");
			scanf("%f", &t);
			if (t==1 || t==-1)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][5]=t;
		
		//ToDo: gripper after Move Here
		good=0;
		while(!good)
		{
			printf("Gripper after move(Open = 1, Close = -1): ");
			scanf("%f", &t);
			if (t!=0)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][6]=t;
		
		good=0;
		while(!good)
		{
			printf("Enter speed of move(0-9): ");
			scanf("%f", &t);
			if (-.01<=t && t<=9.01)
			{
				good=1;
			}
			else
				printf("Value out of range. Try again!\n");	
		}
		commands[x][7]=t;
		printf("\n");

		// Start DAMAGE CHECK

		// convert angles to radians
		x2=toRadians(commands[x][1]);
		x3=toRadians(commands[x][2]);
		x4=toRadians(commands[x][3]);

		z5=-10.0*(30.0 + 25.0*sin(x2) + 16.0*sin(x2 + x3));
		z6=-300.0 - 250.0*sin(x2) - 160.0*sin(x2 + x3) - 179.0*sin(x2 + x3 + x4);

		printf("%.1f %.1f\n", z5, z6);

		// test angles 
		if (z5>-75 || z6>-20)
		{
			x--;
			printf("ERROR: Position would result in damage to robot!\n"); 
			printf("Please re-enter position.\n");
          
			// reset fields to zero
			commands[x][1]=0;
			commands[x][2]=0;
			commands[x][3]=0;
		}
	}
	numPostions += i;
}

double toRadians(double degang)
{
	return (degang*3.14159)/180.0;
}

double toDegrees(double rad)
{
	return (rad*180)/3.14159;
}

/**********************************************************************/ 
void MOD_ANGLES () 
/* This subroutine should be written so that the user can modify any given 
          position in memory. Again verify that changes meet the proper ranges, 
          and also run you crash protection. This means you will also need to 
          use CALC_TRANSFORMS() in this routine. */
{
	char buf[2];
	int row, col, good, i;
	float t;
	float oldValue;
	double x2, x3, x4;  
	double z5, z6;      
	
	printf("Enter position number to modify: ");
    scanf("%d", &row);

	// display position
	i = row - 1;
	printf("\nCOL/   1        2        3        4        5       6     7     8  \n"); 
	printf("ROW    BA       SH       EL       WP       WR      GB    GA    SP  \n"); 			
	printf("%1d %8.1f %8.1f %8.1f %8.1f %8.1f %5.0f %5.0f %5.0f", i+1, 
		commands[i][0], commands[i][1], commands[i][2], commands[i][3],
		commands[i][4], commands[i][5], commands[i][6], commands[i][7]);
		printf("\n\n");

	printf("Enter column number to modify: ");
	scanf("%d", &col);
	
	oldValue = commands[row-1][col-1];
	good=0;
	while(!good) 
	{
		printf("Enter new data for position %d column %d: ", row, col);
		scanf("%f", &t);
		if (col == 1) 
		{
			if (-150<=t && t<=150)
			{
				good=1;
			}
			else
				printf("Incorrect value. BA (-150 to 150)\n");	
		}
		
				//Todo: col == 2
	
		if (col == 2) 
		{
			if (-30<=t && t<=105)
			{
				good=1;
			}
			else
				printf("Incorrect value. BA (-30 to 105)\n");	
		}
				//Todo: col == 3
	
		if (col == 3) 
		{
			if (-105<=t && t<=2)
			{
				good=1;
			}
			else
				printf("Incorrect value. BA (-105 to 2)\n");	
		}
				//Todo: col == 4
	
		if (col == 4) 
		{
			if (-90<=t && t<=90)
			{
				good=1;
			}
			else
				printf("Incorrect value. BA (-90 to 90)\n");	
		}
				//Todo: col == 5
		
		if (col == 5) 
		{
			if (-180<=t && t<=180)
			{
				good=1;
			}
			else
				printf("Incorrect value. BA (-180 to 180)\n");	
		}
						
		
		if (col == 6)
		{
			if (t==1 || t==-1)
			{
				good=1;
			}
			else
				printf("Incorrect value. GB (Open = 1, Close = -1)\n");	
		}
		
		//Todo: col == 7
		if (col == 7)
		{
			if (t==1 || t==-1)
			{
				good=1;
			}
			else
				printf("Incorrect value. GB (Open = 1, Close = -1)\n");	
		}
		
		if (col == 8)
		{
			if (-.01<=t && t<=9.01)
			{
				good=1;
			}
			else
				printf("Incorrect value. SP (0 to 9)\n");	
		}
		commands[row - 1][col - 1] = t;

		// Start DAMAGE CHECK

		// convert angles to radians
	
		x2=toRadians(commands[row-1][1]);
		x3=toRadians(commands[row-1][2]);
		x4=toRadians(commands[row-1][3]);

		z5=-10.0*(30.0 + 25.0*sin(x2) + 16.0*sin(x2 + x3));
		z6=-300.0 - 250.0*sin(x2) - 160.0*sin(x2 + x3) - 179.0*sin(x2 + x3 + x4);
				
		// test angles 
		if (z5>-75 || z6>-20)
		{
			printf("\nERROR: Position would result in damage to robot!\n"); 
			printf("Please re-enter position.\n");
          
			// reset to old value
			commands[row-1][col-1] = oldValue;
			good = 0;
		}
	}
	printf("Change Successful!\n\n");
	printf("COL/   1        2        3        4        5       6     7     8  \n"); 
	printf("ROW    BA       SH       EL       WP       WR      GB    GA    SP  \n"); 			
	printf("%1d %8.1f %8.1f %8.1f %8.1f %8.1f %5.0f %5.0f %5.0f", i+1, 
		commands[i][0], commands[i][1], commands[i][2], commands[i][3],
		commands[i][4], commands[i][5], commands[i][6], commands[i][7]);
		printf("\n\n");
	
	// hold screen display
	scanf("%c", &buf);
	printf ("\nPress <enter> to continue... "); 
	getchar();

}


/**********************************************************************/ 
void CALC_TRANSFORMS(double x1, double x2, double x3, double x4, double x5) 
/* This is where you'll calculate T05 and T06. Simply rewrite your 
           Mathematica results for those two matricies into this routine. 
          You will need to pass the five joint angles (in radians) to this 
           subroutine. The results should be stored in an array of some type 
           so that the calculation is not performed over and over. */
{ 
    double l1, l2, l3, lt;

    //Link lengths       
    l1 = 300;
	l2 = 250;
	l3= 160;
	lt = -179;
           
           	// T05 transformation calculations 
	
	t5[1][1] = C(x1)*C(x5)*S(x2+x3+x4)-S(x1)*S(x5);
	t5[1][2] = -C(x5)*S(x1)-C(x1)*S(x2+x3+x4)*S(x5);
	t5[1][3] = -C(x1)*C(x2+x3+x4);
	t5[1][4]= C(x1)*(l2*C(x2)+l3*C(x2+x3));
	t5[2][1] = -C(x5)*S(x1)*S(x2+x3+x4)-C(x1)*S(x5);
	t5[2][2]= -C(x1)*C(x5)+S(x1)*S(x2+x3+x4)*S(x5);
	t5[2][3] = C(x2+x3+x4)*S(x1);
	t5[2][4]= -(l2*C(x2)+l3*C(x2+x3))*S(x1);
	t5[3][1]= -C(x2+x3+x4)*C(x5);
	t5[3][2] = C(x2+x3+x4)*S(x5);
	t5[3][3] = -S(x2+x3+x4);
	t5[3][4] = l1+l2*S(x2)+l3*S(x2+x3);
	t5[4][1]=0;
	t5[4][2]=0;
	t5[4][3]=0;
	t5[4][4]=1;
           
    // T06 transformation calculations 
	
	t6[1][1] = C(x1)*C(x5)*S(x2+x3+x4)-S(x1)*S(x5);
	t6[1][2] = -C(x5)*S(x1)+C(x1)*S(x2+x3+x4)*S(x5);
	t6[1][3] = C(x1)*C(x2+x3+x4);
	t6[1][4] = C(x1)*(l2*C(x2)+l3*C(x2+x3)-lt*C(x2+x3+x4));
	t6[2][1] = -C(x5)*S(x1)*S(x2+x3+x4)-C(x1)*S(x5);
	t6[2][2]= C(x1)*C(x5)-S(x1)*S(x2+x3+x4)*S(x5);
	t6[2][3]= -C(x2+x3+x4)*S(x1);
	t6[2][4]= -(l2*C(x2)+l3*C(x2+x3)-lt*C(x2+x3+x4))*S(x1);
	t6[3][1]= -C(x2+x3+x4)*C(x5);
	t6[3][2] = -C(x2+x3+x4)*S(x5);
	t6[3][3] = S(x2+x3+x4);
	t6[3][4] = l1+l2*S(x2)+l3*S(x2+x3)-lt*S(x2+x3+x4);
	t6[4][1]=0;
	t6[4][2]=0;
	t6[4][3]=0;
	t6[4][4]=1;
/****           
system("clear"); 
printf ("\n\n\nNo positions exist.\n");***/
 
printf ("\nPress <enter> to continue... "); 
getchar(); 
}


/*********************************************************************/ 
/* This simply prints out you stored positions to the screen. */

void PRINT_POS () 
{
	char buf[2];
	int p = numPostions;
    int i; 
	printf("COL/   1        2        3        4        5       6     7     8  \n"); 
	printf("ROW    BA       SH       EL       WP       WR      GB    GA    SP  \n"); 			
 
	for (i=0; i<p; i++) 
	{
		printf("%1d %8.1f %8.1f %8.1f %8.1f %8.1f %5.0f %5.0f %5.0f", i+1, 
			commands[i][0], commands[i][1], commands[i][2], commands[i][3],
				commands[i][4], commands[i][5], commands[i][6], commands[i][7]);
			printf("\n");
	}
	scanf("%c", &buf);
   	printf ("\nPress <enter> to continue... "); 
	getchar();
}

/*********************************************************************/ 
void printTransforms() 
/* This subroutine prints the T05 and T06 matricies for a given position 
           to the screen. Try to make it as neat as possible since this is likely 
           to be the screen I observe most of the time. */
{ 
// task: display position and orientation of tool frame w.r.t. base frame

	char buf[2];
	int pos, p;

	double x1, x2, x3, x4, x5;
	/*double r11, r12, r13, r14;
	double r21, r22, r23, r24;
	double r31, r32, r33, r34;
	double r41, r42, r43, r44;
	double q11, q12, q13, q14;
	double q21, q22, q23, q24;
	double q31, q32, q33, q34;
	double q41, q42, q43, q44;*/
	double l1, l2, l3, lt;

	l1 = 300;
	l2 = 250;
	l3= 160;
	lt = -179;

	printf("Enter an existing position number: ");
	scanf("%d", &pos);
    	
	p = pos - 1;
    int cont = 1;
    while(cont)
    {
	// convert degrees to radians for transformation calculations	
	x1=toRadians(commands[p][0]);
	x2=toRadians(commands[p][1]);
	x3=toRadians(commands[p][2]);
	x4=toRadians(commands[p][3]);
	x5=toRadians(commands[p][4]);

	printf("test of x1 = %.1f x2 = %.1f x3 = %.1f x4 = %.1f x5 = %.1f\n", x1, x2, x3, x4, x5);

	// BEGINNING-OF-TO5 ***********************************************

     CALC_TRANSFORMS(x1,x2,x3,x4,x5);
	// T05 transformation calculations 
	
/**	r11 = C(x1)*C(x5)*S(x2+x3+x4)-S(x1)*S(x5);
	r12 = -C(x5)*S(x1)-C(x1)*S(x2+x3+x4)*S(x5);
	r13 = -C(x1)*C(x2+x3+x4);
	r14 = C(x1)*(l2*C(x2)+l3*C(x2+x3));
	r21 = -C(x5)*S(x1)*S(x2+x3+x4)-C(x1)*S(x5);
	r22= -C(x1)*C(x5)+S(x1)*S(x2+x3+x4)*S(x5);
	r23 = C(x2+x3+x4)*S(x1);
	r24= -(l2*C(x2)+l3*C(x2+x3))*S(x1);
	r31= -C(x2+x3+x4)*C(x5);
	r32 = C(x2+x3+x4)*S(x5);
	r33 = -S(x2+x3+x4);
	r34 = l1+l2*S(x2)+l3*S(x2+x3);
	r41=0;
	r42=0;
	r43=0;
	r44=1;
	**/
	
	
	// T05 display transformation matrix 
	printf("      _                                               _ \n");
	printf("     | %10.4f %10.4f %10.4f %12.4f   |\n", t5[1][1], t5[1][2], t5[1][3], t5[1][4]); 
	printf("0    | %10.4f %10.4f %10.4f %12.4f   |\n", t5[2][1], t5[2][2], t5[2][3], t5[2][4]);
	printf(" T = | %10.4f %10.4f %10.4f %12.4f   |\n", t5[3][1], t5[3][2], t5[3][3], t5[3][4]);
	printf("5    |_%10.4f %10.4f %10.4f %12.4f  _|\n", t5[4][1], t5[4][2], t5[4][3], t5[4][4]);
	printf("\n");

	// END-OF-T05 *****************************************************


	// BEGINNING-OF-TO6 ***********************************************

	// T06 transformation calculations 
	
/**	q11 = C(x1)*C(x5)*S(x2+x3+x4)-S(x1)*S(x5);
	q12 = -C(x5)*S(x1)+C(x1)*S(x2+x3+x4)*S(x5);
	q13 = C(x1)*C(x2+x3+x4);
	q14 = C(x1)*(l2*C(x2)+l3*C(x2+x3)-lt*C(x2+x3+x4));
	q21 = -C(x5)*S(x1)*S(x2+x3+x4)-C(x1)*S(x5);
	q22=  C(x1)*C(x5)-S(x1)*S(x2+x3+x4)*S(x5);
	q23 = -C(x2+x3+x4)*S(x1);
	q24= -(l2*C(x2)+l3*C(x2+x3)-lt*C(x2+x3+x4))*S(x1);
	q31= -C(x2+x3+x4)*C(x5);
	q32 = -C(x2+x3+x4)*S(x5);
	q33 = S(x2+x3+x4);
	q34 = l1+l2*S(x2)+l3*S(x2+x3)-lt*S(x2+x3+x4);
	q41=0;
	q42=0;
	q43=0;
	q44=1;
  **/ 	
	// T06 display transmation matrix 
	printf("      _                                               _ \n");
	printf("     | %10.4f %10.4f %10.4f %12.4f   |\n", t6[1][1], t6[1][2],t6[1][3],t6[1][4]); 
	printf("0    | %10.4f %10.4f %10.4f %12.4f   |\n", t6[2][1], t6[2][2], t6[2][3], t6[2][4]);
	printf(" T = | %10.4f %10.4f %10.4f %12.4f   |\n", t6[3][1], t6[3][2], t6[3][3], t6[3][4]);
	printf("6    |_%10.4f %10.4f %10.4f %12.4f  _|\n", t6[4][1], t6[4][2], t6[4][3], t6[4][4]);
	printf("\n");

	// END-OF-T06 *****************************************************
 
    if(t5[3][4]<0 || t6[3][4]<0)
    {
	    printf(" Position will result in damage to robot - please modify\n");
      
     }
     else
     {
         cont=0;
     }
 

    }
	// hold screen display
	scanf("%c", &buf);
	printf ("\nPress <enter> to continue... "); 
	getchar();
}

/*************************************************************************/ 
void MOVE_ROBOT () 
/* This subroutine allows the user to view the current positions, then to 
           select any number of positions, in any given order, to move the robot 
           to. 
*/ 
{

//RVM1_MoveJoint( float WaistAngle,
//			         float ShoulderAngle,
//			         float ElbowAngle,
//			         float WristPitchAngle,
//			         float WristRotationAngle )
	int moves[100];
	float next[5];
	int t, x, y, max=0;
   
	printf("How many positions for the loop: ");
	scanf("%d", &max);
	PRINT_POS();
    
	for(x=0; x<max; x++)
	{
		printf("Enter position #: ");
		scanf("%d", &t);
		moves[x]=t-1;
	}
    RVM1_Origin();
	for(x=0; x<max; x++)
	{
		t=moves[x];
		if (commands[t][5]>0)
		{
			RVM1_GripOpen();
		}
		else
		{
			RVM1_GripClose();
		}
		RVM1_Speed((int)commands[t][7]);
		
		for(y=0; y<5; y++)
		{
			next[y]=commands[t][y]-curr[y];
		}

		RVM1_MoveJoint(next[0], next[1], next[2], next[3], next[4]);

		for(y=0; y<5; y++)
		{
			curr[y]=commands[t][y];
		}
		if (commands[t][6]>0)
		{
			RVM1_GripOpen();
		}
		else
		{
			RVM1_GripClose();
		}

		getchar();
		
	}


}

/************************************************************************/ 
void GRASP_OBJECT() 
{
	getchar();

}
/************************************************************************/
double C(double degrees)
{      
      return cos(degrees);
       
}
/***********************************************************************/
double S(double degrees)
{
      
       return sin(degrees);
}
