/** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
/**                                                             **/
/**  DO NOT MODIFY THIS FILE!  DOING SO MAY IMPACT YOUR         **/
/**    PROGRAM'S ABILITY TO COMMUNICATE WITH THE RV-M1 ROBOT    **/
/**    CONTROLLER!!                                             **/
/**                                                             **/
/** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/

/* RobotSerialInterface.c */

#include <windows.h>
#include <stdio.h>
#include <winbase.h>
#include <string.h>

#include	"RobotSerialInterface.h"


HANDLE PortHandle;

FILE*	RVM1_LogFile;
static unsigned int    RVM1_Options;


/** ----------------------------------------------------------- **/
unsigned int RVM1_InitializeComm( unsigned int InitializeOptions,
									char * LogFile)
{
	char LogFileName[255];
	char commport[10];  /*name of serial port*/

	extern HANDLE PortHandle;
	
	DCB CommSettings;
	COMMTIMEOUTS timeout;

    printf("\n\nRobot Commands Will Be Sent To:\n");
    
    if( InitializeOptions & RVM1_WriteToFile )
	{
		if( LogFile == "" )
		{
			strcpy( LogFileName, "C:\\_RVM1_Logs\\RVM1_CommLog.txt" );
		}
		else
		{
			strcpy( LogFileName, LogFile );
		}

		if( (RVM1_LogFile = fopen( LogFileName, "w+" )) == NULL )
		{
			printf( "The file %s was not opened\n", LogFileName );
		}
		else
		{
			printf( "\tLog File %s\n", LogFileName );
			RVM1_Options = RVM1_Options | RVM1_WriteToFile;
		}
	}

	if( InitializeOptions & RVM1_WriteToScreen )
	{
        printf("\tConsole Screen\n");
			RVM1_Options = RVM1_Options | RVM1_WriteToScreen;
    }

	if( InitializeOptions & RVM1_WriteToComm )
	{
        printf("\tRV-M1 Robot Controller via COM port\n");
        /**  mode com1:4800,n,8,1,p  **/

		/* Creating a handle to the port */
		sprintf(commport,"COM1");
		PortHandle = CreateFile(commport, GENERIC_WRITE | GENERIC_READ,
						0, NULL, OPEN_EXISTING, 0, NULL);
		
		if( SetupComm( PortHandle, 2048, 2048) == 0 )
		{
			printf("\n\nBuffer Sizes NOT Updated\n		");
		}

		/* Settings for commport */
		GetCommState(PortHandle, &CommSettings);
		CommSettings.DCBlength = sizeof(DCB);
		CommSettings.BaudRate = CBR_4800;
		CommSettings.fBinary = TRUE;
		CommSettings.fParity = FALSE;
		CommSettings.fOutxCtsFlow = TRUE;
		CommSettings.fOutxDsrFlow = TRUE;
		CommSettings.fDtrControl = DTR_CONTROL_HANDSHAKE;
		CommSettings.fDsrSensitivity = FALSE;
		CommSettings.fTXContinueOnXoff = FALSE;
		CommSettings.fOutX = TRUE;
		CommSettings.fInX = FALSE;
		CommSettings.fErrorChar = FALSE;
		CommSettings.fNull = FALSE;
		CommSettings.fRtsControl = RTS_CONTROL_HANDSHAKE;
		CommSettings.fAbortOnError = FALSE;
		CommSettings.ByteSize = 8;
		CommSettings.Parity = NOPARITY;
		CommSettings.StopBits = ONESTOPBIT;
		SetCommState(PortHandle, &CommSettings);

		/* Settings timeout values */
		GetCommTimeouts(PortHandle,&timeout);
		timeout.ReadIntervalTimeout = 500;
		timeout.ReadTotalTimeoutMultiplier = 500;
		timeout.ReadTotalTimeoutConstant = 500;
		timeout.WriteTotalTimeoutMultiplier = 1000;
		timeout.WriteTotalTimeoutConstant = 1000;
		SetCommTimeouts(PortHandle, &timeout);

		RVM1_Options = RVM1_Options | RVM1_WriteToComm;
	}
    printf("\n\n");
	return RVM1_Options;
}


/** ----------------------------------------------------------- **/
void RVM1_CleanUp( void )
{
	if( RVM1_Options & RVM1_WriteToFile )
	{
        fclose(RVM1_LogFile);
        printf("RVM1 Log File closed ....\r\r");
    }

    fflush( stdout );
    fflush( stdin );
    fflush( stderr );
}


/** ----------------------------------------------------------- **/
void RVM1_FlushBuffers( void )
{
	if( RVM1_Options & RVM1_WriteToFile )
	{
        fflush( RVM1_LogFile );
    }

	if( RVM1_Options & RVM1_WriteToComm )
	{
        FlushFileBuffers( PortHandle );
    }
}






/**                                                             **/
/**  Robot Specific Commands                                    **/
/** ----------------------------------------------------------- **/
void RVM1_GripClose( void )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[25];

	ReturnStatus = sprintf(RobotCommand, "GC\r\n");

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

	if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);
    }
}


/** ----------------------------------------------------------- **/
void RVM1_GripOpen( void )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[25];

	ReturnStatus = sprintf(RobotCommand, "GO\r\n");

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

	if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);
    }
}


/** ----------------------------------------------------------- **/
void RVM1_MoveJoint( float WaistAngle,
			         float ShoulderAngle,
			         float ElbowAngle,
			         float WristPitchAngle,
			         float WristRotationAngle )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[50];

	ReturnStatus = sprintf(RobotCommand, 
                    "MJ %3.1f, %3.1f, %3.1f, %3.1f, %3.1f\r\n",
			            WaistAngle,
			            ShoulderAngle,
			            ElbowAngle,
			            WristPitchAngle,
			            WristRotationAngle);

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

    if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);
    }
}


/** ----------------------------------------------------------- **/
void RVM1_Nest( void )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[25];

	ReturnStatus = sprintf(RobotCommand, "NT\r\n");

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

	if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);	
	}
}


/** ----------------------------------------------------------- **/
void RVM1_Origin( void )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[25];

	ReturnStatus = sprintf(RobotCommand, "OG\r\n");

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

	if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);
    }
}


/** ----------------------------------------------------------- **/
void RVM1_Reset( void )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[25];

	ReturnStatus = sprintf(RobotCommand, "RS\r\n");

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

	if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);
    }
}


/** ----------------------------------------------------------- **/
void RVM1_Speed( int MoveSpeed )
{
	int		BytesWritten;
	int		ReturnStatus;
	char	RobotCommand[25];

    if( MoveSpeed > 9 )
    {
        MoveSpeed = 9;
    }

    if( MoveSpeed < 0 )
    {
        MoveSpeed = 0;
    }

	ReturnStatus = sprintf(RobotCommand, "SP %i\r\n", MoveSpeed);

	if( RVM1_Options & RVM1_WriteToScreen )
	{
        printf("\n\n  **RV-M1 DEBUG ->  %s",RobotCommand);
    }

	if( RVM1_Options & RVM1_WriteToFile )
	{
        fprintf(RVM1_LogFile,"%s",RobotCommand);
    }
    
	if( RVM1_Options & RVM1_WriteToComm )
	{
        /**  Send Command to RV-M1  **/
	    /** --------------------------------------------------  **/
	    WriteFile(PortHandle,
			      RobotCommand,
			      strlen(RobotCommand),
			      &BytesWritten,
			      NULL);
    }
}



