//originally in main
//

/* string getDirName(const string& prefix, const Domain& dom) */ 
/* { */
/*   stringstream filenameSS; */
/*   filenameSS << "out/" << prefix << "_" << dom.Ni << "/"; */
/*   return filenameSS.str(); */
/* } */    
/* #include <config5cpp/Configuration.h> */
/* using namespace config4cpp; */
void getState(Vector3d &u, Matrix3d &F, double &rho, double& S, Setting&);

int readFromFile()
{
	/* Config cfg; */
	/* Matrix3d F_t; */
	/* Vector3d u_t; */
	/* double S_t, rho_t, e_t; */
	/* try */
	/* { */
	/* 	cfg.readFile("./Configurations/testCases.cfg"); */
	/* } */
	/* catch(const FileIOException &fioex) */
	/* { */
	/* 	std::cerr << "I/O error while reading file." << std::endl; */
	/* 	return(EXIT_FAILURE); */
	/* } */
	/* catch(const ParseException &pex) */
	/* { */
	/* 	std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine() */
	/* 		<< " - " << pex.getError() << std::endl; */
	/* 	return(EXIT_FAILURE); */
	/* } */
	/* try */
	/* { */
	/* 	cfg.lookupValue("slChoice", slLimChoice); */
	/* 	cfg.lookupValue("base", baseDir); */
	/* 	Setting &settingState = cfg.lookup("testCases.testCase1.left"); */
	/* 	cout << "Length " << settingState.getLength() << endl; */
	/* 	int temp = settingState["rho"].getLength(); */

	/* 	//get leftstate */
	/* 	getState(u_t, F_t, rho_t, S_t, settingState); */
	/* 	e_t = internalEnergy(F_t, S_t); */	
	/* 	ElasticState ElasticStateL(u_t, F_t, rho_t, S_t, e_t); */
	/* 	cout << ElasticStateL << endl; */	
	/* 	settingState.lookupValue("rho", rho_t); */
	/* 	cout << rho_t << endl; */
	/* 	cout << "slChoice" <<slLimChoice << endl; */
	/* } */
	/* catch(const SettingNotFoundException &nfex) */
	/* { */
	/* 	cerr << "No 'name' setting in configuration file." << endl; */
	/* } */
}

// get settings from configuration file
void getState(Vector3d &u, Matrix3d &F, double &rho, double& S, Setting& setting)
{
	setting.lookupValue("rho", rho);
	setting.lookupValue("S", S);
	int max = setting["u"].getLength();
	double temp = 0; 
	for (unsigned int i = 0; i < max; i++)
	{
		u(i) = setting["u"][i];
	}
	int count = 0;
	for (unsigned int i = 0; i < max; i++)
	{
		for (unsigned int j = 0; j < max; j++)
		{
			count = i*max + j;
			F(i,j) = setting["F"][count];
		}
	}
}
