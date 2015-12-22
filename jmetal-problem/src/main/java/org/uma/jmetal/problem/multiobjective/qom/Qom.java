package org.uma.jmetal.problem.multiobjective.qom;

import org.uma.jmetal.problem.ConstrainedProblem;
import org.uma.jmetal.problem.impl.AbstractDoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.solutionattribute.impl.OverallConstraintViolation;

import java.awt.*;
import java.awt.List;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Class representing problem Qom
 * Simple hydrological model for the calculation of losses and excesses
 */
public class Qom extends AbstractDoubleProblem implements ConstrainedProblem<DoubleSolution> {

  /**
   * Constructor.
   * Creates a default instance of the Qom problem.
   * is a Hydrological Model Rainfall-Runoff
   * Gustavo R. Zavala
   * Antonio J. Nebro.
   * @param solutionType The solution type must "Real" or "BinaryReal".
   */

  // selected objetive functions
  protected  String  file_;
  protected  String[] selectedOF_;
  // amount evaluations
  protected  int  numberOfEval_;
  // maximum evaluations
  protected  int test_;

  // Basin file
  protected  String description_;
  protected  String shapeName_; // Fast, Triangle, Trapeze, Rectangle, Oval, Slow
  protected  double shape_; // coefficient of bsin shape, 0.5 to 2.5
  protected  double area_; // km2
  protected  double pervious_; // 0.0 to 1.0
  protected  double CT_; // concentration time
  protected String dateIni_ ; // "yyyy/MM/dd hh:mm:00 tt"
  protected String dateFin_ ; // "yyyy/MM/dd hh:mm:00 tt"
  protected  double Δt_; // hours decimal
  protected  int regIni_;
  protected  int regFin_;
  protected Boolean percolationIsBaseFlow_; // True or False
  protected  double Qoi_; // initial flow in impervious area (m3/s)
  protected  double Qop_; // initial flow in pervious area (m3/s)
  protected  double Qog_; // initial flow in ground water (m3/s)
  protected  double Ro_; // humedad antecedente volume (mm)
  protected  double So_; // humedad antecedente volume (mm)
  protected  double Rmax_; // mm
  protected  double Smax_; // mm
  protected  double PC_; // mm/h
  protected  double Ki_; // 1/s
  protected  double Kp_; // 1/s
  protected  double Kg_; // 1/s

  protected  int varRmax_=0;
  protected  int varSmax_=1;
  protected  int varPC_=2;
  protected  int varKi_=3;
  protected  int varKp_=4;
  protected  int varKg_=5;
  protected  int varCT_=6;
  protected  int varShape_=7;

  //conditions wet initial
  double Ri_1_;
  double Si_1_;
  double Di_1_;
  // potential volumes
  double Dp_;
  // total volumes
  double sumVi_;
  double sumVp_;
  double sumVg_;
  // moisture balance
  double Input_;
  double Evapotraspiration_;
  double SuperficialStore_;
  double SoilStore_;
  double StoreOrOutput_Deep_ ;
  double OutputSuperficial_ ;
  double MassBalance_;
  double MassBalanceError_;
  double Qpo_;// maximum streamflow for hydrograph observaded
  double Qpc_;// maximum streamflow for hydrograph calculated
  double VHo_;// volume for hydrograph observaded
  double VHc_;// volume for hydrograph calculated
  double maxO_;// maximum observed steamflow

  String textOF_="";
  String optOF_="";
  double optRmax_=0;
  double optSmax_=0;
  double optPc_=0;
  double optKi_=0;
  double optKp_=0;
  double optKg_=0;
  double optCT_=0;
  double optShape_=0;

  // state variable for read VAR file
  int maxLine_, line_ ;
  // read variables in VAR jMetal file
  double[][] var_ ;
  // objetive function derivaded for validation test
  double[][] fd_;

  //  rainfall (constant)
  protected double [] Pconst_ ;
  public double getcconst_(int reg) {
    return Pconst_[reg] ;
  } // get Pconst

  //  rainfall accumulated (constant)
  protected double [] Paccum_ ;
  public double getPaccum_(int reg) {
    return Paccum_[reg] ;
  } // get Paccum

  //  rainfall
  protected double [] P_ ;
  public double getP_(int reg) {
    return P_[reg] ;
  } // get P

  //  potential evapotraspiration (constant)
  protected double [] PEVTconst_ ;
  public double getPEVTconst_(int reg) {
    return PEVTconst_[reg] ;
  } // get PEVTconst

  //  potential evapotraspiration
  protected double [] PEVT_ ;
  public double getPEVT_(int reg) {
    return PEVT_[reg] ;
  } // get PEVT

  //  total rainfall
  protected double [] EVT_ ;
  public double getEVT_(int reg) {
    return EVT_[reg] ;
  } // get EVT

  //  observed flow
  protected double [] Qo_ ;
  public double getQo_(int reg) {
    return Qo_[reg] ;
  } // get Qo

  //  calculated flow
  protected double [] Qc_ ;
  public double getQc_(int reg) {
    return Qc_[reg] ;
  } // get Qc

  //  observed depth
  protected double [] Ho_ ;
  public double getHo_(int reg) {
    return Ho_[reg] ;
  } // get Ho

  //  calculated depth
  protected double [] Hc_ ;
  public double getHc_(int reg) {
    return Hc_[reg] ;
  } // get Hc

  //  superficial storage
  protected double [] R_ ;
  public double getR_(int reg) {
    return R_[reg] ;
  } // get R

  //  soil storage
  protected double [] S_ ;
  public double getS_(int reg) {
    return S_[reg] ;
  } // get S

  //  deep storage
  protected double [] D_ ;
  public double getD_(int reg) {
    return D_[reg] ;
  } // get D

  //  impervious volume for time i
  protected double [] Vi_ ;
  public double getVi_(int reg) {
    return Vi_[reg] ;
  } // get Vi

  //  pervious volume for time i
  protected double [] Vp_ ;
  public double getVp_(int reg) {
    return Vp_[reg] ;
  } // get Vp

  //  ground water volume for time i
  protected double [] Vg_ ;
  public double getVg_(int reg) {
    return Vg_[reg] ;
  } // get Vg

  //  impervious flux for time i
  protected double [] Qi_ ;
  public double getQi_(int reg) {
    return Qi_[reg] ;
  } // get Qi

  //  pervious flux for time i
  protected double [] Qp_ ;
  public double getQp_(int reg) {
    return Qp_[reg] ;
  } // get Qp

  //  ground water flux for time i
  protected double [] Qg_ ;
  public double getQg_(int reg) {
    return Qg_[reg] ;
  } // get Qg

  //  impervious volume for time i
  protected double [] Vis_ ;
  public double getVis_(int reg) {
    return Vis_[reg] ;
  } // get Vis

  //  pervious volume for time i
  protected double [] Vps_ ;
  public double getVps_(int reg) {
    return Vps_[reg] ;
  } // get Vps

  //  ground water volume for time i
  protected double [] Vgs_ ;
  public double getVgs_(int reg) {
    return Vgs_[reg] ;
  } // get Vgs

  //  impervious flux for time i
  protected double [] Qis_ ;
  public double getQis_(int reg) {
    return Qis_[reg] ;
  } // get Qis

  //  pervious flux for time i
  protected double [] Qps_ ;
  public double getQps_(int reg) {
    return Qps_[reg] ;
  } // get Qps

  //  ground water flux for time i
  protected double [] Qgs_ ;
  public double getQgs_(int reg) {
    return Qgs_[reg] ;
  } // get Qgs

  //  volume of hydrograph observated
  protected double [] Vo_ ;
  public double getVo_(int reg) {
    return Vo_[reg] ;
  } // get Vo

  //  volume of hydrograph calculated
  protected double [] Vc_ ;
  public double getVc_(int reg) {
    return Vc_[reg] ;
  } // get Vc

  public OverallConstraintViolation<DoubleSolution> overallConstraintViolationDegree ;

  public Qom() {
    line_=0;
    overallConstraintViolationDegree = new OverallConstraintViolation<DoubleSolution>() ;
    Qom_Initialize();
  }

  public void Qom_Initialize() {

    // name of class and problem
    setName("Qom");

    System.out.println("Qom a simple model of losses and excesses");
    System.out.println("Qom-Clark a rainfall-runoff model");
    // read definition problem
    Qom_ReadProblems("Qom.txt");
    // dataset read files
    System.out.println("Problem file: " + file_);
    Qom_ReadFileBasin(file_ + ".basin");

    //shape_=ShapeCoefficient(shapeName_);
    //Fill lower and upper limits

    Qom_ReadLimitsVariable(file_ + ".var");
    Qom_ReadFileRain(file_ + ".rain");
    Qom_ReadFilePEVT(file_ + ".pevt");
    Qom_ReadFileQo(file_ + ".Qo");
    // Qom_ReadFileHo(file_ + ".Ho");
    // problem data print
    System.out.println("Optimization multi-objective: ");
    System.out.println("  Number of objective function: " + getNumberOfObjectives());
    System.out.println("  " + textOF_);
    System.out.println("  Number of Variables: " + getNumberOfVariables());
    System.out.println("  Number of Constraints: " + getNumberOfConstraints());

    System.out.println("Algorithm configuration: ");

    optCT_= CT_;
    optShape_=shape_;
    numberOfEval_= 0;

  } // end InitializeEBEs

  public void evaluate(DoubleSolution solution) {

    int hi=0;
    double [] fx = new double[getNumberOfObjectives()] ; // functions

    if (test_ == -1 || test_ == 1)
    {
      if(test_ == -1 && line_ < maxLine_) // variables from VAR
      {
        System.out.println("Evaluation: " +  line_) ;
        for(int j=0;j<getNumberOfVariables();j++){
          solution.setVariableValue(j, var_[line_][j]);
        }
        if(Ro_ > solution.getVariableValue(varRmax_))
          Ro_=solution.getVariableValue(varRmax_);
        if(So_ > solution.getVariableValue(varSmax_))
          So_=solution.getVariableValue(varSmax_);
      }
      else if(test_==1)
      {
        // variables form file *.basin
        solution.setVariableValue(varRmax_, Rmax_);
        solution.setVariableValue(varSmax_, Smax_);
        solution.setVariableValue(varPC_, PC_);
        solution.setVariableValue(varKi_, Ki_);
        solution.setVariableValue(varKp_, Kp_);
        solution.setVariableValue(varKg_, Kg_);
        if(getNumberOfVariables()>=7)
          solution.setVariableValue(varCT_, CT_);
        if(getNumberOfVariables()>=8)
          solution.setVariableValue(varShape_, shape_);

        if(Ro_ > Rmax_ || So_ > Smax_ ){
          System.out.println("Error: antecedent moisture conditions are greater than the");
          System.out.println("allowable maximum, check the values ​​of Ro, So and limit variables");
          System.exit(1);
        }
      }
    }
    else{
      if(Ro_ > solution.getVariableValue(varRmax_))
        solution.setVariableValue(varRmax_, Ro_);
      if(So_ > solution.getVariableValue(varSmax_))
        solution.setVariableValue(varSmax_, So_);
    }

    Rmax_=solution.getVariableValue(varRmax_);
    Smax_=solution.getVariableValue(varSmax_);
    PC_=solution.getVariableValue(varPC_); //
    Ki_=solution.getVariableValue(varKi_); //
    Kp_=solution.getVariableValue(varKp_); //
    Kg_=solution.getVariableValue(varKg_); //
    if(getNumberOfVariables()>=7) CT_=solution.getVariableValue(varCT_);
    if(getNumberOfVariables()>=8) shape_=solution.getVariableValue(varShape_);

    // Method simple losses, a hydrological problem
    Qom_Calculus();

    // START OBJETIVES FUNCTION
    for (int i = 0; i < getNumberOfObjectives(); i++) {
      if (selectedOF_[i].equals("NSE")) {
        fx[i] = FunctionNSE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("NSE-V")) {
        fx[i] = FunctionNSE(Vo_, Vc_);
      } else if (selectedOF_[i].equals("NSEO")) {
        fx[i] = FunctionNSEO(Qo_, Qc_);
      } else if(selectedOF_[i].equals("NSEO-V")){
        fx[i] = FunctionNSEO(Vo_, Vc_);
      } else if (selectedOF_[i].equals("R2")) {
        fx[i] = FunctionCDR2(Qo_, Qc_);
      } else if (selectedOF_[i].equals("PBIAS")) {
        fx[i] = FunctionPBIAS(Qo_, Qc_);
      } else if (selectedOF_[i].equals("RQV")) {
        fx[i] = FunctionRQV(Qo_, Qc_, Vo_, Vc_);
      } else if(selectedOF_[i].equals("ERQV")) {
        fx[i] = FunctionERQV(Qo_, Qc_, Vo_, Vc_);
      } else if (selectedOF_[i].equals("Module")) {
        fx[i] = FunctionModule(Qo_, Qc_);
      } else if (selectedOF_[i].equals("SE")) {
        fx[i] = FunctionSE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("SIE")) {
        fx[i] = FunctionSIE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("RSE")) {
        fx[i] = FunctionRSE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("MSE")) {
        fx[i] = FunctionMSE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("M4E")) {
        fx[i] = FunctionM4E(Qo_, Qc_);
      } else if (selectedOF_[i].equals("MSDE")) {
        fx[i] = FunctionMSDE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("MLSE")) {
        fx[i] = FunctionMLSE(Qo_, Qc_);
      } else if (selectedOF_[i].equals("EMQQ")) {
        fx[i] = FunctionEMQQ(Qo_, Qc_);
      } else if(selectedOF_[i].equals("MBE")) {
        double a = Paccum_[P_.length-1];
        double b = (sumVg_+sumVi_+sumVp_);
        fx[i] = (1-Math.min(a,b)/Math.max(a,b))*100;
      } else if (selectedOF_[i].equals("REMB")) {
        // Relative Error of the Mass Balance
        fx[i] = Math.abs(MassBalance_/Input_)*100;
      } else {
        System.out.println("Error: the objective function " + i + " is not recognized");
        System.exit(1);
      }

      solution.setObjective(i, fx[i]);
    }
    // END OBJETIVES FUNCTION

    VHo_=SumVolumeOfHydrograph(Vo_);
    VHc_=SumVolumeOfHydrograph(Vc_);

    numberOfEval_++;
  } // evaluate

  /**
   * Evaluates the constraint overhead of a solution
   * @param solution The solution
   * @throws JMetalException
   */
  public void evaluateConstraints(DoubleSolution solution) {

    double [] constraint = new double[this.getNumberOfConstraints()];
    double [] fx = new double[getNumberOfObjectives()] ; // functions
    for(int i=0;i<getNumberOfObjectives();i++)
    {
      fx[i]=solution.getObjective(i);
    }

    // constraint No 1: Mass Balance
    double absoluteVolumeErr = MassBalance_;
    double observedVolume = Input_;
    double relativeVolumeErr =absoluteVolumeErr / observedVolume;
    constraint[0] = 0.30 -  Math.abs(relativeVolumeErr);// < 0.25 -> 0.25%

    // constraint No 2: Nash-Sutclife Efficiency
    double NSE = FunctionNSE(Qo_, Qc_);
    constraint[1] = NSE - 0.70; // < at 0.85 of NSE

    // constraint No 3: relative error for volume derived hydrograph
    // < at 25.0->0.25%
    // x [Hm3] * 1000.0 [mm/m] / Area[km2]
    observedVolume = (VHo_*1000.0/area_);
    if(percolationIsBaseFlow_ && true)
    {
      absoluteVolumeErr= (StoreOrOutput_Deep_ + OutputSuperficial_) - observedVolume;
    }
    else{
      absoluteVolumeErr = OutputSuperficial_ - observedVolume;
    }
    relativeVolumeErr = absoluteVolumeErr / observedVolume;
    constraint[2] = 0.30 - Math.abs(relativeVolumeErr); // < at 0.25->25%

    double total = 0.0;
    int number = 0;
    String optOF ="";
    for (int i = 0; i < this.getNumberOfConstraints(); i++)
      if (constraint[i]<0.0){
        total+=constraint[i];
        number++;
      }
      else{
        // objective functions
        // decision variables
        optRmax_=solution.getVariableValue(varRmax_);
        optSmax_=solution.getVariableValue(varSmax_);
        optPc_=solution.getVariableValue(varPC_);
        optKi_=solution.getVariableValue(varKi_);
        optKp_=solution.getVariableValue(varKp_);
        optKg_=solution.getVariableValue(varKg_);
        if (getNumberOfVariables()>=7)
          optCT_=solution.getVariableValue(varCT_);
        if (getNumberOfVariables()>=8)
          optShape_=solution.getVariableValue(varShape_);

        //System.out.println("Evaluation: " + numberOfEval_) ;
        optOF ="";
        for(int j=0;j<getNumberOfObjectives();j++) {
          optOF += "f" + (j+1) + ": " + solution.getObjective(j)+ " ";
        }

        if((numberOfEval_ % 1000) == 0)
        {
          System.out.println("Evaluation: " + numberOfEval_) ;
          System.out.println(optOF);
          // System.out.println("  Volumes hydrograph. O: " + VHo_ + " - C: " + VHc_ );
          System.out.println(" Rmax: " + optRmax_); //x[0]
          System.out.println(" Smax: " + optSmax_); //x[1]
          System.out.println(" PC: " + optPc_); //x[2]
          System.out.println(" Ki: " + optKi_); //x[3]
          System.out.println(" Kp: " + optKp_); //x[4]
          System.out.println(" Kg: " + optKg_); //x[5]
          System.out.println(" CT: " + optCT_); //x[6]
          System.out.println(" Shape : " + optShape_); //x[7]
        }
      }

    //solution.setOverallConstraintViolation(total);
    //solution.setNumberOfViolatedConstraint(number);
    overallConstraintViolationDegree.setAttribute(solution, total);


    if(test_ == 1 || test_ == -1)
    {
      optOF ="";
      for(int j=0;j<getNumberOfObjectives();j++) {
        optOF += "f" + (j+1) + ": " + solution.getObjective(j)+ " ";
      }
      System.out.println(optOF); //x[0]
      System.out.println(" Rmax: " + optRmax_); //x[0]
      System.out.println(" Smax: " + optSmax_); //x[1]
      System.out.println(" PC: " + optPc_); //x[2]
      System.out.println(" Ki: " + optKi_); //x[3]
      System.out.println(" Kp: " + optKp_); //x[4]
      System.out.println(" Kg: " + optKg_); //x[5]
      System.out.println(" CT: " + optCT_); //x[6]
      System.out.println(" Shape : " + optShape_); //x[7]

      if (test_ == 1)
      {
        Qom_PrintArchSolutions();
        if(test_ == 1) System.exit(0);
      }
      else if (test_ == -1)
      {
        fd_[line_][0]=fx[0];
        fd_[line_][1]=fx[1];
        if(line_>= maxLine_-1) {
          Qom_save_OF_validation();
          System.exit(0);
        }
        line_++;
      }
    }
    else {
      if(numberOfEval_==test_)
        Toolkit.getDefaultToolkit().beep();
    }

  } // evaluateConstraints

  public void Qom_Calculus() throws JMetalException {

    // clear variables

    Qom_InitializeProblemVariables();

    // initials condition
    // previous storage
    Ri_1_ = Ro_;
    double z = Math.exp(-0.0001 * Δt_ * Ro_ / (Rmax_ - Ro_) * pervious_);
    Si_1_ = Math.max(Smax_ * (1.0 - z), So_);
    Di_1_ = Math.min(PC_ * Δt_ * (1.0 - z), Si_1_);

    for (int i = 0; i < P_.length; i++)
    {
      // evaporation and traspiration potetencial
      // and interception include in surface, vegetation.

      if (Ri_1_ >= PEVT_[i])
      {
        LossesReservoir_Excess(i);
        ExcessRain(i);
        SoilStorage(i);
        LossesPercolation(i);
      }
      else
      {
        LossesReservoir_Deficit(i);

        if (P_[i] > PEVT_[i])
        {
          RechargesRain(i);
          ExcessRain(i);
          SoilStorage(i);
          LossesPercolation(i);
        }
        else
        {
          LossesRain(i);
          LossesSoil(i);
          LossesPercolation(i);
        }
      }

      // separation of volumes
      SeparationVolumes(i);

      Ri_1_ = R_[i];
      Si_1_ = S_[i];
      Di_1_ = D_[i];
    }

    Runoff();

    MoistureBalance();
  }

  public void Qom_InitializeProblemVariables() throws JMetalException {

    // problem variables
    P_ = new double[regFin_ - regIni_ + 1];
    Paccum_ = new double[regFin_ - regIni_ + 1];
    PEVT_ = new double[regFin_ - regIni_ + 1];
    EVT_ = new double[regFin_ - regIni_ + 1];
    Qc_ = new double[regFin_ - regIni_ + 1];
    Hc_ = new double[regFin_ - regIni_ + 1];
    R_ = new double[regFin_ - regIni_ + 1];
    S_ = new double[regFin_ - regIni_ + 1];
    D_ = new double[regFin_ - regIni_ + 1];
    Vi_ = new double[regFin_ - regIni_ + 1];
    Vp_ = new double[regFin_ - regIni_ + 1];
    Vg_ = new double[regFin_ - regIni_ + 1];
    Qi_ = new double[regFin_ - regIni_ + 1];
    Qp_ = new double[regFin_ - regIni_ + 1];
    Qg_ = new double[regFin_ - regIni_ + 1];
    Vis_ = new double[regFin_ - regIni_ + 1];
    Vps_ = new double[regFin_ - regIni_ + 1];
    Vgs_ = new double[regFin_ - regIni_ + 1];
    Qis_ = new double[regFin_ - regIni_ + 1];
    Qps_ = new double[regFin_ - regIni_ + 1];
    Qgs_ = new double[regFin_ - regIni_ + 1];
    Vo_ = new double[regFin_ - regIni_];
    Vc_ = new double[regFin_ - regIni_];


    double Plast=0.0;
    for (int i = 0; i < Pconst_.length; i++)
    {
      P_[i]= Pconst_[i];
      Paccum_[i]=Plast+P_[i];
      PEVT_[i] = PEVTconst_[i];

      Plast = Paccum_[i];
    }
  }

  public void LossesReservoir_Excess(int i)
  {
    R_[i] = Ri_1_ - PEVT_[i];
    EVT_[i] = PEVT_[i];
    PEVT_[i] = 0.0;
  }

  public void LossesReservoir_Deficit(int i)
  {
    PEVT_[i] -= Ri_1_;
    EVT_[i] = Ri_1_;
    R_[i] = 0.0;
  }

  public void LossesSoil(int i)
  {
    // deficit soil stogae
    // line relaation
    // S_[i] = Math.max(Si_1_ * (1.0 - PEVT_[i] / Smax_), 0.0);
    double z2 = Math.exp(-0.05 * Si_1_ / (Smax_ - Si_1_));
    S_[i] = Math.max(Si_1_ - PEVT_[i] * z2, 0.0);
    double ΔS = Si_1_ - S_[i];
    EVT_[i] +=  ΔS;
    PEVT_[i] -= ΔS;
    // capillary rise
    double CR = Di_1_ * 0.20; //0.10
    if (CR > ΔS && (S_[i] + CR <=Smax_ ))
    {
      S_[i] = S_[i] + CR;
      Di_1_ *= 0.80; //0.90
    }
    Dp_ = Di_1_;

  }

  public void LossesRain(int i)
  {
    // rain is less than the evapotranspiration
    PEVT_[i] -= P_[i];
    EVT_[i] += P_[i];
    P_[i] = 0.0;
  }

  public void RechargesRain(int i)
  {
    // is P >= PEPT
    EVT_[i] += PEVT_[i];
    P_[i] -= PEVT_[i];
    PEVT_[i] = 0.0;
  }

  public void ExcessRain(int i)
  {
    if (P_[i] > Rmax_ - R_[i])
    {
      // rain is greater than the rest of the storage
      P_[i] = P_[i] - Rmax_ + R_[i];
      R_[i] = Rmax_;
    }
    else
    {
      // rain is less than the storage
      R_[i] +=  P_[i];
      P_[i] = 0.0;
    }
  }

  public void SoilStorage(int i)
  {
    double z = Math.exp(-0.0001 * Δt_ * R_[i] / (Rmax_ - R_[i]) * pervious_);
    double Sp = Smax_ * (1.0 - z);
    S_[i]=Math.max(Si_1_, Sp);
    double ΔS = Math.max(S_[i] - Si_1_, 0.0);

    if (P_[i] - ΔS * 0.40 >= 0.0) //0.40
    {
      P_[i] -= ΔS * 0.40;
    }
    else
    {
      ΔS = P_[i];
      P_[i] = 0.0;
      S_[i] = Math.min(Si_1_ + ΔS, Smax_);
    }

    R_[i] = Math.max(R_[i] - ΔS, 0.0);
    Dp_ = PC_ * Δt_ * (1 - z); //
  }

  public void LossesPercolation(int i)
  {
    if (S_[i] >= Dp_)
    {
      D_[i] = Dp_;
      S_[i] = S_[i] - D_[i];
    }
    else
    {
      D_[i] = 0.0;
    }
  }

  public void SeparationVolumes(int i)
  {
    // water volume in impervious area
    Vi_[i] = P_[i] * (1.0-pervious_);
    // water volume in pervious area no saturated
    Vp_[i] = P_[i] * pervious_;
    // deep groudwater volume
    Vg_[i] = D_[i];
  }

  public void Runoff() throws JMetalException {
    // surface spatial distribution coefficient applied to surface runoff volumes
    Vis_=SyntheticHistogramArea(shape_, CT_, Δt_, Vi_);
    Vps_=SyntheticHistogramArea(shape_, CT_, Δt_, Vp_);
    if (percolationIsBaseFlow_)
    {
      Vgs_=SyntheticHistogramArea(shape_, CT_, Δt_, Vg_);
    }

    // overlandwater flux
    Qis_=SimpleLinearReservoir(area_, Vis_, Ki_, Δt_, Qoi_);

    //groundwater layer
    Qps_=SimpleLinearReservoir(area_, Vps_, Kp_, Δt_, Qop_);

    // deep water flow for the groundwater
    if (percolationIsBaseFlow_)
    {
      Qgs_=SimpleLinearReservoir(area_, Vgs_, Kg_, Δt_, Qog_);
    }
    // total flux spread calculated
    for (int i = 0; i < Qis_.length; i++)
    {
      Qc_[i] = Qis_[i] + Qps_[i] + Qgs_[i];
    }

    Vo_=VolumeOfHydrograph(Qo_);
    Vc_=VolumeOfHydrograph(Qc_);
  }

  public void MoistureBalance()
  {
    Input_=0.0;
    Evapotraspiration_=0.0;
    sumVi_=0.0;
    sumVp_=0.0;
    sumVg_=0.0;

    int numReg=Pconst_.length;
    for(int i=0;i<numReg;i++)
    {
      Input_+=Pconst_[i];
      Evapotraspiration_ += EVT_[i];
      sumVi_+=Vi_[i];
      sumVp_+=Vp_[i];
      sumVg_+=Vg_[i];
    }

    SuperficialStore_ = R_[numReg-1];
    SoilStore_ = S_[numReg-1];
    StoreOrOutput_Deep_ = sumVg_;
    OutputSuperficial_ = sumVi_ + sumVp_;
    // estimated - observed
    MassBalance_ = -So_ -Ro_ -Input_ + SuperficialStore_ + SoilStore_ + StoreOrOutput_Deep_ + OutputSuperficial_ + Evapotraspiration_;
    // mass balance errors
    MassBalanceError_ = MassBalance_ / Input_ * 100;
  }

  public double ShapeCoefficient(String shape)  throws JMetalException
  {
    /// <summary>
    /// shape coefficient
    /// </summary>
    double n=1.0;
    try {
      // coeficient shape number
      // n is watershed response, if n is less than 1 quick, if n is greater than 1 is slow
      shape = shapeName_.toUpperCase();
      if (shape.equals("FAST")) {
        n = 0.25;
      } else if (shape.equals("TRIANGLE")) {
        n = 0.5;
      } else if (shape.equals("TRAPEZE")) {
        n = 0.75;
      } else if (shape.equals("RECTANGLE")) {
        n = 1.0;
      } else if (shape.equals("OVAL")) {
        n = 1.25;
      } else if (shape.equals("DIAMOND")) {
        n = 1.50;
      } else if (shape.equals("SLOW")) {
        n = 1.75;
      } else {
        System.out.println("Error: Shape Area Coefficient not recognized ( " + shape + " )");
        System.exit(1);
      };
    }
    catch (Exception ex) {
      System.out.println("Shape Area Coefficient - Error: " + ex.getMessage());
    }

    return n;
  }

  public double[] SyntheticHistogramArea(double n, double CT, double Δt, double[] V)  throws JMetalException
  {
    /// <summary>
    /// shape coefficient number
    /// </summary>
    /// <summary>
    /// CT concentration time (hour)
    /// </summary>
    /// <summary>
    ///  Δt step time (hour)
    // / </summary>
    /// <summary>
    ///  V excess volume (mm)
    // / </summary>

    /// <summary>
    /// cumulative area
    /// </summary>
    double[] Ac;
    /// <summary>
    /// area between two consecutive cumulative areas
    /// </summary>
    double[] f;
    /// <summary>
    /// number of subbasin subdivision for for distribution
    /// </summary>
    int No;

    double[] Vs = new double[V.length];
    try
    {
      // coeficient
      //double a = Math.pow(0.5, (1.0 - n)); // de prueba
      double a = Math.pow(0.5, n); // original

      // number division of cuenca
      No = (int)(CT / Δt);
      if (No < 1)
      {
        No = 1;
      }

      // determining the cumulative areas
      Ac = new double[No];
      double[] time = new double[No];
      double dt = 1.0/No;
      for(int i=0;i<No;i++)
      {
        time[i]=(i+1)*dt;
      }
      // t: is this unitary concentration time
      // 1/No: is this step of time
      for(int i = 0; i<No;i++)
      {
        if (time[i] <= 0.5)
        {
          Ac[i] = a * Math.pow(time[i], n);
        }
        else
        {
          Ac[i] = 1.0 - a * Math.pow((1.0 - time[i]), n);
        }
      }

      // area between two consecutive cumulative areas
      f = new double[No];
      f[0] = Ac[0];
      for (int i = 1; i < No; i++)
      {
        f[i] = Ac[i] - Ac[i - 1];
      }

      // runoff volume spead in function of Synthetic Histogram Area
      // limit down
      double fc = 0.0f;

      // runoff volume spread
      for (int t = 0; t < V.length; t++)
      {
        for (int i = 0; i <= t; i++)
        {
          if (t - i < No) fc = f[t - i]; else fc = 0.0;
          // ratio for spread and sum volume runoff imprevius
          Vs[t] += V[i] * fc;
          // Console.WriteLine("{0}({1},{2})={3}", t, i, t-i, Vs[t]);
        }
      }
    }
    catch (Exception ex) {
      System.out.println("SynteticHistogramAreas - Error: " + ex.getMessage());
    }

    return Vs;
  }

  public double[] SimpleLinearReservoir(double Area, double[] Vs, double K, double Δt, double Qo)  throws JMetalException
  {
    /// spread flow undergoes a buffer pool of equivalent
    /// to a simple linear function.
    /// overland water using pervious and impervious areas and groundwater
    /// Area: area basin (km2)
    /// Vs: volume of the streamflow in (mm/h)
    /// Ki, Kp and Kg: average time (hour) of emptying the reservoir related to Vs
    /// Δt: step time (hours)
    /// Qo: initial flux (m3/s)
    /// Qs: streamflow spread (m3/s)

    double[] Qs = new double [Vs.length];

    try
    {
      Qs[0]=Qo;

      //conversion mm/h to m3/s for area catchement
      // streamflow [m3/s] = V[mm/h] * 1[h/3600s] * 1000000[m2/km2] * Area[Km2] * 1[m/1000mm]
      Double coefConvert = area_ * 1.0 / 3.6 / Δt_;

      if (K > 0.0)
      {
        Qs[0] = Qs[0] + (Vs[0] * coefConvert * (1.0f - Math.pow(Math.E, (-Δt / K))));

        for (int i = 1; i < Vs.length; i++)
        {
          Qs[i] = Qs[i - 1] * Math.pow(Math.E, (-Δt / K)) + Vs[i] * coefConvert * (1.0f - Math.pow(Math.E, (-Δt / K)));
        }
      }
    }
    catch (Exception ex)
    {
      System.out.println("Simple linear Reservoir - Error: " + ex.getMessage());
    }

    return Qs;
  }

  public double FunctionNSE(double[] O, double[] E)
  {
    // function Efficiency Nash-Sutcliffe
    // O[k] : k-th observed data value of streamflow
    // E[k] : k-th estimated or forecasted value of streamflow
    int k = O.length;
    double Om = 0.0;

    //mean of the observed data
    for(int i=0; i<k;i++)
    {
      Om += O[i];
    }
    Om = Om/k;

    double SSRes = 0.0;
    double SSTot = 0.0;
    double NSE = 0.0;
    for (int i = 0; i < k; i++)
    {
      //Sum of Squares of Residuals, also called the residual sum of squares
      SSRes += Math.pow((E[i] - O[i]), 2.0);
      //Total Sum of Squares (proportional to the sample variance)
      SSTot += Math.pow((O[i] - Om), 2.0);
    }
    NSE = 1.0 - SSRes / SSTot;
    return NSE;
  }

  public double FunctionNSEO(double[] O, double[] E)
  {
    // function Efficiency Nash-Sutcliffe
    // O[k] : k-th observed data value of streamflow
    // E[k] : k-th forecasted or estimated or calculated of value of streamflow (discharge)

    double NSE = 0.0;
    double NSEO = 1.0;
    NSE = FunctionNSE(O, E);
    NSEO = 1.0 - NSE;

    return NSEO;
  }

  public double FunctionPBIAS(double[] O, double[] F)
  {
    // function percent bias (PBIAS)
    // measures the average tendency of the simulated data to be larger
    // or smaller than their observed counterparts
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted or estimated or calculated of value of streamflow (discharge)
    double PBIAS=0.0;
    double SOMF=0.0;
    double SO=0.0;

    // register numbers
    int k = O.length;
    for(int i=0; i<k;i++)
    {
      //if(test_>1) {
      //    if (maxO_ > ) F[i] *= 0.85;
      //}

      SOMF += (F[i]  - O[i]);
      SO += O[i];
    }

    PBIAS = SOMF / SO * 100.0;

    return Math.abs(PBIAS);
  }

  public double FunctionSE(double[] O, double[] F)
  {
    // function Square Error
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted or estimated or calculated of value of streamflow (discharge)
    int k = O.length;
    double SE=0.0;
    for(int i=0; i<k;i++)
    {
      SE += Math.pow((F[i] - O[i]), 2.0 );
    }
    return SE;
  }

  public double FunctionSIE(double[] O, double[] F)
  {
    // function Square Inverse Error
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    int k = O.length;
    double SIE=0.0;
    for(int i=0; i<k;i++)
    {
      if(O[i]!=0.0 && F[i]!=0.0 && (1.0/F[i] - 1.0/O[i]) > 0.0)
      {
        SIE += Math.pow((1.0/F[i] - 1.0/O[i]), 2.0 );
      }
    }
    return SIE;
  }

  public double FunctionModule(double[] O, double[] F)
  {
    // function Module
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    int k = O.length;
    double M=0.0;
    for(int i=0; i<k;i++)
    {
      M += Math.abs(F[i] - O[i]);
    }
    return M;
  }

  public double FunctionRSE(double[] O, double[] F)
  {
    // function Relative Square Error
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    double RSE=0.0;
    // register numbers
    int k = O.length;
    for(int i=0; i<k;i++)
    {
      RSE += Math.pow((F[i] - O[i])/O[i], 2.0 );
    }
    return RSE;
  }

  public double FunctionMSE(double[] O, double[] F)
  {
    // function Mean Error Square-root
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    double MSE=0.0;
    // register numbers
    int k = O.length;
    for(int i=0; i<k;i++)
    {
      MSE += Math.sqrt(Math.abs(F[i] - O[i]));
    }
    MSE  = MSE / k;
    return MSE;
  }

  public double FunctionRQV(double[] Qo, double[] Qf, double[]Vo, double[]Vf)
  {
    // function Ratio Streamflow Volume
    // Qo[k] : k-th observed data value of streamflow
    // Qf[k] : k-th forecasted value of streamflow
    // Vo[k] : l-th observed data value of streamflow
    // Vf[k] : l-th forecasted value of streamflow
    double RQV=0.0;
    // register numbers
    int k = Qo.length;
    int l = Vo.length;

    double Qom = 0;
    double RQoc = 0;
    double RQom = 0.0;
    double Vom = 0.0;
    double RVoc = 0;
    double RVom = 0.0;

    for(int i=0; i<k;i++)
    {
      Qom += Qo[i];
    }
    Qom = Qom/k;

    for(int j=0; j<l;j++)
    {
      Vom += Vo[j];
    }
    Vom = Vom/l;

    for (int i = 0; i < k; i++)
    {
      RQoc += Math.pow((Qf[i] - Qo[i]), 2.0);
      RQom += Math.pow((Qo[i] - Qom), 2.0);
    }

    for (int j = 0; j < l; j++)
    {
      RVoc += Math.pow((Vf[j] - Vo[j]), 2.0);
      RVom += Math.pow((Vo[j] - Vom), 2.0);
    }
    double a = RQoc / RQom;
    double b = RVoc / RVom;
    RQV =  1- b / a; //

    return RQV;
  }

  public double FunctionCDR2(double[] O, double[] E)
  {
    // function Coefficient of Determination
    // denoted R2 and pronounced R squared
    // O[k] : k-th observed data value of streamflow
    // E[k] : k-th estimated or forecasted value of streamflow
    int N = O.length; // total No of data set
    double SO = 0.0; // Sum of observed data
    double SE = 0.0; // sum of values estimated
    double SOE = 0.0;
    double SO2 = 0.0;
    double SE2 = 0.0;
    double R=0.0;
    double R2=0.0;

    //mean of the observed data
    for(int i=0; i<N;i++)
    {
      SO += O[i];
      SE += E[i];
      SOE += O[i]*E[i];
      SO2 += O[i]*O[i];
      SE2 += E[i]*E[i];
    }
    // Pearson’s correlation coefficient
    R = (N*SOE-SO*SE)/(Math.sqrt((N*SO2-Math.pow(SO,2.0))*(N*SE2-Math.pow(SE,2.0))));
    // coefficient of determination
    R2= R*R;

    return R2;
  }

  public double FunctionERQV(double[] Qo, double[] Qf, double[]Vo, double[]Vf)
  {
    // function enhanced Ratio Streamflow and Volume
    // Qo[k] : k-th observed data value of streamflow
    // Qf[k] : k-th forecasted value of streamflow
    // Vo[k] : l-th observed data value of streamflow
    // Vf[k] : l-th forecasted value of streamflow
    double RQV=0.0;

    // register numbers
    int k = Qo.length;
    int l = Vo.length;

    double Qom = 0;
    double RQoc = 0;
    double RQom = 0.0;
    double Vom = 0.0;
    double RVoc = 0;
    double RVom = 0.0;

    for(int i=0; i<k;i++)
    {
      Qom += Qo[i];
    }
    Qom = Qom/k;

    for(int j=0; j<l;j++)
    {
      Vom += Vo[j];
    }
    Vom = Vom/l;

    for (int i = 0; i < k; i++)
    {
      RQoc += Math.pow((Qf[i] - Qo[i]), 2.0);
      RQom += Math.pow((Qo[i] - Qom), 2.0);
    }

    for (int j = 0; j < l; j++)
    {
      RVoc += Math.pow((Vf[j] - Vo[j]), 2.0);
      RVom += Math.pow((Vo[j] - Vom), 2.0);
    }
    //RQV = Math.sqrt(RVoc / Vo.length) * RQom / RQoc;

    double m = 1.0;
    int p=(int)(RVoc * m);
    while (p < 1)
    {
      m *= 10.0;
      p=(int)(RVoc * m);
    }
    RQV = Math.sqrt(RVoc * m / Vo.length) * RQom / RQoc * 100.0;

    return RQV;
  }

  public double FunctionMLSE(double[] O, double[] F)
  {
    // function mean logarithmic squared error (MLSE)
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    double MSLE=0.0;
    int k=O.length;
    for(int i=0; i<k;i++)
    {
      if(O[i] != 0.0 && F[i] != 0.0)
      {
        MSLE += Math.pow(Math.log(F[i])-Math.log(O[i]), 2.0);
      }
    }
    MSLE = MSLE/k;

    return MSLE;
  }

  public double FunctionMSDE(double[] O, double[] F)
  {
    // function mean squared derivative error (MSDE)
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    double MSDE=0.0;
    int k=O.length;
    for(int i=1; i<k;i++)
    {
      MSDE += Math.pow((F[i]-F[i-1])-(O[i]-O[i-1]),2.0);
    }
    MSDE = MSDE/(k-1);

    return MSDE;
  }

  public double FunctionM4E(double[] O, double[] F)
  {
    // mean fourth-power error (M4E) (de Vos and Rientjes, 2007, 2008)
    // O[k] : k-th observed data value of streamflow
    // F[k] : k-th forecasted value of streamflow
    double M4E=0.0;
    int k=O.length;
    for(int i=1; i<k;i++)
    {
      M4E += Math.pow((F[i]-O[i]), 4.0);
    }
    M4E = M4E/k;

    return M4E;
  }

  public double FunctionEMQQ(double[] O, double[] F)
  {
    double EMQQ=0.0;
    double minQo=1.0E208;
    double minQf=1.0E208;
    double maxQo=-1.0E208;
    double maxQf=-1.0E208;

    for(int n=0;n< Qo_.length;n++)
    {
      if(minQo>=Qo_[n])
      {
        minQo=Qo_[n];
        minQf=Qc_[n];
      }

      if(maxQo<=Qo_[n])
      {
        maxQo=Qo_[n];
        maxQf=Qc_[n];
      }
    }
    double EQmin = Math.sqrt(Math.pow((minQf - minQo), 2.0));
    double EQmax = Math.sqrt(Math.pow((maxQf - maxQo), 2.0));
    EMQQ = EQmin+EQmax;

    return EMQQ;
  }

  public double[] VolumeOfHydrograph(double[] f)
  {
    // f: function in m3/s
    // for volume convertUnit to Hm3 = X[m3/s]*Δt[h]*3600[s/h]*1[Hm3/1000000] = Δt * 0.36 / 100.0
    // for flows convertUnit (m3/s) = 1.0
    double convertUnit = Δt_*0.36/100.0;
    double[] fm = new double[f.length-1];
    for (int j = 1; j < f.length; j++)
    {
      fm[j-1] = (f[j - 1] + f[j]) / 2.0 * convertUnit;
    }
    return fm;
  }

  public double SumVolumeOfHydrograph(double[] f)
  {
    double area=0.0;
    for (int j = 0; j < f.length; j++)
    {
      area+=f[j];
    }

    return area;
  }

  public void maxStreamflows()
  {
    Qpo_=-1.0e208;
    Qpc_=-1.0e208;
    for (int j = 1; j < Qo_.length; j++)
    {
      if(Qpo_ < Qo_[j]) Qpo_ = Qo_[j];
      if(Qpc_ < Qc_[j]) Qpc_ = Qc_[j];
    }
  }

  public void Qom_ReadProblems(String fileName) throws JMetalException{

    char ch;
    int lineN = 0;
    String line = "";
    String var1 = "";
    String var2 = "";
    String fileLoad = fileName;
    int i=0, j=0;
    selectedOF_= new String[2];

    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Head read
      lineN++;
      line=input.nextLine();
      lineN++;
      line=input.nextLine();
      lineN++;
      line=input.nextLine();

      // file read
      line=input.nextLine();
      j=0;
      for(i=line.length()-1;i>=0;i--){
        ch=line.charAt(i);
        if(ch == ' '){
          j=i+1;
          break;
        }
      }
      test_ = Integer.valueOf(line.substring(j));

      // number of decision Variables
      int w=0;
      for(i=j-2;i>=0;i--){
        ch=line.charAt(i);
        if(ch == ' '){
          w=i+1;
          break;
        }
      }
      setNumberOfVariables(Integer.valueOf(line.substring(w, j-1)));

      // number of objective function
      int v=0;
      for(i=w-2;i>=0;i--){
        ch=line.charAt(i);
        if(ch == ' '){
          v=i+1;
          break;
        }
      }
      setNumberOfObjectives(Integer.valueOf(line.substring(v,w-1)));

      // number of constraint function
      setNumberOfConstraints(3);

      selectedOF_= new String[getNumberOfObjectives()];
      int m = 0;
      for(j=getNumberOfObjectives()-1;j>-1;j--) {
        for (i = v - 2; i >= 0; i--) {
          ch = line.charAt(i);
          if (ch == ' ') {
            m = i + 1;
            break;
          }
        }
        selectedOF_[j] = line.substring(m, v - 1);
        v=m;
        m=0;
      }

      textOF_="";
      for(i=0;i<getNumberOfObjectives();i++) {
        textOF_ += "f" + (i+1) + ": " + selectedOF_[i]+ " ";
      }

      file_ = line.substring(0, v-1);

      // close the file
      input.close();

      // for validation at VAR jMetal file
      if(test_==-1)
      {
        // count amount lines
        fileLoad = "ReferenceFrontAndVariableStudy-jMetal.txt";
        file = new java.io.File(fileLoad);
        Scanner s = new Scanner(file);
        maxLine_=-1; // head
        do{
          String linea = s.nextLine();
          maxLine_++;
        }while (s.hasNextLine());
        // close the file
        input.close();

        // read variables in VAR jMetal file
        fd_ = new double [maxLine_][getNumberOfObjectives()];
        var_ = new double [maxLine_][getNumberOfVariables()];
        file = new java.io.File(fileLoad);
        s = new Scanner(file);
        String linea = s.nextLine(); // head
        i=0;
        do{
          linea = s.nextLine();
          String[] campos = linea.split(" ");
          for(String campo : campos){
            int ind=0;
            for(j=getNumberOfObjectives();j<(getNumberOfObjectives()+getNumberOfVariables());j++)
            {
              var_[i][ind]=Double.valueOf(campos[j]);
              ind++;
            }
          }
          i++;
        }while (s.hasNextLine());
        // close the file
        input.close();
      }
    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileLoad + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(0);
    }
  }

  public void Qom_ReadFileBasin(String fileName) throws JMetalException{

    int i, j=0, lineN=0;
    char ch;
    String txt = "";

    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Read data from file

      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();

      lineN++;
      // description
      txt=input.nextLine();
      for(i=0;i<txt.length();i++){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      description_ = txt.substring(j);

      lineN++;
      // shape
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      //shapeName_ = txt.substring(j);
      txt = txt.substring(j);
      shape_= Double.valueOf(txt);

      lineN++;
      // number of bar groups
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      area_=Double.valueOf(txt);

      lineN++;
      // pervious
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      pervious_=Double.valueOf(txt);

      lineN++;
      // concentration time
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      CT_=Double.valueOf(txt);

      lineN++;
      // initial date-time
      txt=input.nextLine();
      txt = txt.substring(txt.length()-24);
      dateIni_= txt;

      lineN++;
      // final date-time
      txt=input.nextLine();
      txt = txt.substring(txt.length()-24);
      dateFin_=txt;

      lineN++;
      // time step
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Δt_=Double.valueOf(txt);

      lineN++;
      // initial register
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      regIni_=Integer.valueOf(txt);

      lineN++;
      // final register
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      regFin_=Integer.valueOf(txt);

      lineN++;
      // percolation Is Base Flow ?
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      percolationIsBaseFlow_=Boolean.valueOf(txt);

      lineN++;
      // initial flow in impervious soil
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Qoi_=Double.valueOf(txt);

      lineN++;
      // initial flow in pervious soil
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Qop_=Double.valueOf(txt);

      lineN++;
      // initial downwater
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Qog_=Double.valueOf(txt);

      lineN++;
      // superficial storage initial (reservoir)
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Ro_=Double.valueOf(txt);

      lineN++;
      // soil storage initial (reservoir)
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      So_=Double.valueOf(txt);

      lineN++;
      // superficial storage maximum (reservoir)
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Rmax_=Double.valueOf(txt);

      lineN++;
      // soil storage maximum
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Smax_=Double.valueOf(txt);

      lineN++;
      // ratio ercolation
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      PC_=Double.valueOf(txt);

      lineN++;
      // basin coeficient retention impervious soil
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Ki_=Double.valueOf(txt);

      lineN++;
      // basin coeficient retention pervious soil
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Kp_=Double.valueOf(txt);

      lineN++;
      // basin coeficient retention internal soil
      txt=input.nextLine();
      for(i=txt.length()-1;i>=0;i--){
        ch=txt.charAt(i);
        if(ch == ' '){
          j=i+1;break;}
      }
      txt = txt.substring(j);
      Kg_=Double.valueOf(txt);
      // clse the file
      input.close();
      System.out.println("*.basin file has been read");

    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileName + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(1);
    }
  }

  public void Qom_ReadFileRain(String fileName) throws JMetalException{

    int i, j=0, lineN=0;
    char ch;
    String txt = "";

    double Paccum_1 = 0.0;
    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Read data from file

      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();

      for (i=0;i<regIni_-1;i++){
        lineN++;
        txt=input.nextLine();
      }

      Pconst_ = new double[regFin_ - regIni_ + 1];
      Paccum_ = new double[regFin_ - regIni_ + 1];
      for (i=0;i<regFin_ - regIni_;i++){
        lineN++;
        txt=input.nextLine();
        Pconst_[i]=Double.valueOf(txt);
        Paccum_[i] = Paccum_1 + Pconst_[i];
        Paccum_1 = Pconst_[i];
      }
      input.close();
      System.out.println("*.rain file has been read");
    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileName + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(1);
    }
  }


  public void Qom_ReadFilePEVT(String fileName) throws JMetalException{

    int i, j=0, lineN=0;
    char ch;
    String txt = "";

    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Read data from file

      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();

      for (i=0;i<regIni_-1;i++){
        lineN++;
        txt=input.nextLine();
      }

      PEVTconst_ = new double[regFin_ - regIni_ + 1];
      for (i=0;i<regFin_ - regIni_+ 1;i++){
        lineN++;
        txt=input.nextLine();
        PEVTconst_[i]=Double.valueOf(txt);
      }
      input.close();
      System.out.println("*.pevt file has been read");

    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileName + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(1);
    }
  }

  public void Qom_ReadFileQo(String fileName) throws JMetalException{

    maxO_=0; // maximum observed steamflow
    int i, j=0, lineN=0;
    char ch;
    String txt = "";

    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Read data from file

      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();

      for (i=0;i<regIni_-1;i++){
        lineN++;
        txt=input.nextLine();
      }

      Qo_ = new double[regFin_ - regIni_ + 1];
      for (i=0;i<regFin_ - regIni_+ 1;i++){
        lineN++;
        txt=input.nextLine();
        Qo_[i]=Double.valueOf(txt);
        if (maxO_<Qo_[i]) maxO_=Qo_[i];
      }
      System.out.println("*.Qo file has been read");
      // clse the file
      input.close();
    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileName + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(1);
    }
  }

  public void Qom_ReadFileHo(String fileName) throws JMetalException{

    int i, j=0, lineN=0;
    char ch;
    String txt = "";

    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Read data from file

      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();
      lineN++;
      txt=input.nextLine();

      for (i=0;i<regIni_-1;i++){
        lineN++;
        txt=input.nextLine();
      }

      Ho_ = new double[regFin_ - regIni_ + 1];
      for (i=0;i<regFin_ - regIni_+ 1;i++){
        lineN++;
        txt=input.nextLine();
        Ho_[i]=Double.valueOf(txt);
      }
      System.out.println("*.Ho file has been read");

      // clse the file
      input.close();
    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileName + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(1);
    }
  }

  public void Qom_ReadLimitsVariable(String fileName) throws JMetalException{
    char ch;
    String line = "";
    String var1 = "";
    String var2 = "";
    int i, j=0;
    int lineN=0;

    java.util.List<Double> lowerLimit = new ArrayList<>(getNumberOfVariables()) ;
    java.util.List<Double> upperLimit = new ArrayList<>(getNumberOfVariables()) ;

    try {
      // create a File instance
      java.io.File file = new java.io.File(fileName);
      // create a Scanner for the file
      Scanner input = new Scanner(file);
      // Read data from file

      lineN++;
      line=input.nextLine();
      lineN++;
      line=input.nextLine();
      lineN++;
      line=input.nextLine();

      for(int k=0;k<getNumberOfVariables();k++)
      {
        lineN++;
        line=input.nextLine();
        j=0;
        for(i=line.length()-1;i>=0;i--){
          ch=line.charAt(i);
          if(ch == ' '){
            j=i+1;
            break;
          }
        }
        var2 = line.substring(j);
        upperLimit.add(Double.valueOf(var2));

        int m=0;
        for(i=j-2;i>=0;i--){
          ch=line.charAt(i);
          if(ch == ' '){
            m=i+1;
            break;
          }
        }
        var1 = line.substring(m,j);
        lowerLimit.add(Double.valueOf(var1));
      }

      setLowerLimit(lowerLimit);
      setUpperLimit(upperLimit);

      System.out.println("*.var file has been read");
      // clse the fil e
      input.close();
    }
    catch (Exception ex) {
      System.out.println("ERROR in file: " + fileName + " - (" + ex.getMessage() + ") in line No:" + lineN);
      System.exit(1);
    }
  }

  public void Qom_PrintArchSolutions() throws JMetalException{
    // print solution hydrograph
    try {

      maxStreamflows();

      PrintStream ps = new PrintStream("Optimization.txt");
      ps.printf("Objective Functions"); ps.println();
      ps.printf("Eval: " + numberOfEval_); ps.println();
      ps.printf(textOF_); ps.println();
      ps.printf(optOF_); ps.println();ps.println();
      ps.printf("Rmax (mm): %12.3f", optRmax_); ps.println();
      ps.printf("Smax (mm): %12.3f", optSmax_); ps.println();
      ps.printf("pc (mm/h):   %12.3f", optPc_); ps.println();
      ps.printf("Ki    (h):   %12.3f", optKi_); ps.println();
      ps.printf("Kp    (h):   %12.3f", optKp_); ps.println();
      ps.printf("Kg    (h):   %12.3f", optKg_);
      ps.close();

      ps = new PrintStream("Hydrograph.txt");
      // yyyy/mm/dd hh:mm:ss a.m.
      int year, month, date, hour, minute, second;
      String  tt;
      year =  Integer.valueOf(dateIni_.substring(0,4));
      month =  Integer.valueOf(dateIni_.substring(5,7))-1;
      date =  Integer.valueOf(dateIni_.substring(8,10));
      hour =  Integer.valueOf(dateIni_.substring(11,13));
      minute =  Integer.valueOf(dateIni_.substring(14,16));
      second =  Integer.valueOf(dateIni_.substring(17,19));
      tt =  dateIni_.substring(20,24);
      if(hour==12 && (tt.equals("a.m."))) {
        hour = 0;
      }
      else if(tt.equals("p.m.")) {
        hour+=12;
      }

      Calendar cal = Calendar.getInstance();
      cal.set(year, month, date, hour, minute, second);

      int nDate, nHour, nMinute;
      nDate = 0;
      if((Δt_ % 24.0) == 0) nDate = (int)(Δt_/24.0);
      nHour = 0;
      if((Δt_* 60.0 % 60.0)== 0 && Δt_<24.0 && Δt_>=1.0) nHour=(int)(Δt_* 60.0 / 60.0);
      nMinute = 0;
      if(Δt_< 1.0) nMinute = (int) (60.0 * Δt_);
      SimpleDateFormat sdf = new SimpleDateFormat("yyyy/MM/dd HH:mm");

      //ps.printf("Date Time(h)  Qo (m3/s)   Qc (m3/s)"); ps.println();
      for(int i=0; i< Qo_.length; i++){

        ps.printf(sdf.format(cal.getTime()) + "%10.3f %10.3f %10.3f %10.3f", i*Δt_, Qo_[i], Qc_[i], Qgs_[i]); ps.println();

        cal.add(Calendar.DATE, nDate);
        cal.add(Calendar.HOUR, nHour);
        cal.add(Calendar.MINUTE, nMinute);

      } //Next i

      ps.close();

      ps = new PrintStream("Balabce.txt");
      ps.printf("Mass Balance of Qom"); ps.println();
      ps.printf("Rain (mm):                         %8.2f", Input_); ps.println();
      ps.printf("Initial Superficial Storage (mm):  %8.2f", Ro_); ps.println();
      ps.printf("Initial Soil Storage (mm):         %8.2f", So_); ps.println();
      ps.printf("EVT (mm):                                %8.2f", Evapotraspiration_); ps.println();
      ps.printf("Superficial Storage (mm):                %8.2f", SuperficialStore_); ps.println();
      ps.printf("Soil Storage (mm):                       %8.2f", SoilStore_); ps.println();
      if(percolationIsBaseFlow_)
      {
        ps.printf("Deep Storage (mm):                             -"); ps.println();
        ps.printf("Groundwater Runoff (mm):                 %8.2f", StoreOrOutput_Deep_); ps.println();
      }
      else
      {
        ps.printf("Deep Storage (mm):                       %8.2f", StoreOrOutput_Deep_); ps.println();
        ps.printf("Groundwater Runoff (mm):                       -", 0.0); ps.println();
      }
      ps.printf("Superficial Runoff (mm):                 %8.2f", OutputSuperficial_); ps.println();
      double totalVomumeEstimated = Evapotraspiration_ + SuperficialStore_ + SoilStore_ + StoreOrOutput_Deep_ + OutputSuperficial_;
      ps.printf("Total volume estimated (mm):       %8.2f", totalVomumeEstimated); ps.println();
      ps.printf("Balance (AE mm):                   %8.2f", MassBalance_); ps.println();
      ps.printf("Relative Error (RE)%%:              %8.2f", MassBalanceError_); ps.println();
      ps.printf("-----------------------------------------"); ps.println();
      ps.printf("Hydrograph of Clark"); ps.println();
      ps.printf("Flow peak of Hydrograph"); ps.println();
      ps.printf("Observed (m3/s):                 %8.2f", Qpo_); ps.println();
      ps.printf("Estimated (m3/s):                %8.2f", Qpc_); ps.println();
      double RE=0.0;
      //RE= (1.0-Math.min(Qpo_,Qpc_)/Math.max(Qpo_,Qpc_))*100.0;
      RE= (Qpc_ - Qpo_) / Qpo_ * 100.0;
      ps.printf("Relative error (RE)%%:            %8.2f", RE); ps.println();
      ps.printf("Total Volumes of Hydrograph"); ps.println();
      ps.printf("Observed (Hm3):                  %8.2f", VHo_); ps.println();
      ps.printf("Estimated (Hm3):                 %8.2f", VHc_); ps.println();
      ps.printf("Observed (mm):                   %8.2f", VHo_*1000.0/area_); ps.println();
      ps.printf("Estimated (mm):                  %8.2f", VHc_*1000.0/area_); ps.println();
      //error= (1.0-Math.min(VHc_,VHo_)/Math.max(VHo_,VHc_))*100.0;
      RE= (VHc_-VHo_) / VHo_ * 100.0;
      ps.printf("Relative Error (RE)%%:            %8.2f", RE); ps.println();
      ps.printf("-----------------------------------------"); ps.println();
      ps.printf("Calibration indicators"); ps.println();
      double NSE =FunctionNSE(Qo_, Qc_);
      ps.printf("Nash-Sutcliffe Efficiency (NSE): %8.4f", NSE); ps.println();
      double R2 =FunctionCDR2(Qo_, Qc_);
      ps.printf("Determination coefficient (R2):  %8.4f", R2); ps.println();
      double PBIAS =FunctionPBIAS(Qo_, Qc_);
      ps.printf("Percent bias (PBIAS)          :  %8.4f", PBIAS); ps.println();
      ps.close();

      // volume for hydrograph observaded and hudrograph calculated
      ps = new PrintStream("Volumes.txt");
      //ps.printf("Time(h)  Rain (mm)  Losses(mm) Excess(mm)  R(mm)  S(mm)   D(mm)"); ps.println();
      for(int i=0; i<Pconst_.length; i++){
        ps.printf("%10.3f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", i*Δt_,  Pconst_[i], (Pconst_[i]-P_[i]), P_[i], R_[i], S_[i], D_[i]); ps.println();
      } //Next i
      ps.close();

      // volume for hydrograph observaded and hudrograph calculated
      ps = new PrintStream("Volumes Total of Hydrograph.txt");
      //ps.printf("Time(h)  Vo (m3)   Vc (m3)"); ps.println();
      for(int i=0; i<Vo_.length; i++){
        ps.format("%10.3f %10.6f %10.6f", i*Δt_, Vo_[i], Vc_[i]); ps.println();
      } //Next i
      ps.close();

      System.out.println("Solution files saved");
      Toolkit.getDefaultToolkit().beep();

    }// i
    catch (Exception ex) {
      System.out.println("ERROR in save file Hydrograph.txt: (" + ex.getMessage() + ")");
    }
  }

  public void Qom_save_OF_validation() throws JMetalException{
    // print solution hydrograph
    try {
      PrintStream ps = new PrintStream("Optimization.txt");
      ps = new PrintStream(file_);
      // f1 f2
      for(int i=0; i<maxLine_; i++){
        ps.format(fd_[i][0] + " " + fd_[i][1]); ps.println();
      } //Next i
      ps.close();
      System.out.println("saved function objective derived");
      Toolkit.getDefaultToolkit().beep();
    }
    catch (Exception ex) {
      System.out.println("ERROR in save file SolutionsSort.txt: (" + ex.getMessage() + ")");
    }
  }
}
