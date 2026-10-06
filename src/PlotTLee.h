// Plotting of TLee (included by TLee.cxx, as the mcm_*.h files by master_cov_matrix.cxx): the visualization only, the
// calculations are in TLee.cxx.

#include "PlotTLee_val.h"

// x axis segments of the style in the current pad: an axis with the variable values below the bins of each segment
// (labels if flag_axis_labels), separators between the segments, the segment labels and the test description at the top
// if flag_names
static void draw_gof_segments(const TLeeGoFPlotStyle& style, bool flag_axis_labels, double tick_size, double label_size, bool flag_names)
{
  if( style.segments.empty() ) return;
  gPad->Update();
  double xmin = gPad->GetUxmin(), xmax = gPad->GetUxmax(), ymin = gPad->GetUymin();
  double left = gPad->GetLeftMargin(), right = gPad->GetRightMargin(), bot = gPad->GetBottomMargin(), top = gPad->GetTopMargin();
  auto x_ndc = [&](double x) { return left + (x-xmin)/(xmax-xmin)*(1-left-right); };
  int nseg = style.segments.size();

  for(int iseg=0; iseg<nseg; iseg++) {
	const TLeeGoFAxisSegment& seg = style.segments.at(iseg);

	TGaxis *axis = new TGaxis(seg.index_low, ymin, seg.index_hgh, ymin, seg.value_low, seg.value_hgh, style.segment_divisions, "S");
	axis->SetTickSize(tick_size);
	axis->SetLabelSize(flag_axis_labels ? label_size : 0);
	axis->SetLabelFont(42);
	axis->Draw();

	if( iseg>0 ) {
	  TLine *line_seg = new TLine();
	  line_seg->SetLineColor(kGray+1); line_seg->SetLineWidth(2);
	  line_seg->DrawLineNDC(x_ndc(seg.index_low), bot, x_ndc(seg.index_low), 1-top);
	}

	if( flag_names && seg.label!="" ) {
	  TLatex *latex_seg = new TLatex();
	  latex_seg->SetNDC(); latex_seg->SetTextFont(42);
	  double x_center = x_ndc( 0.5*(seg.index_low + seg.index_hgh + 1) );
	  if( nseg>8 ) { latex_seg->SetTextAngle(90); latex_seg->SetTextAlign(32); latex_seg->SetTextSize(0.03); }
	  else { latex_seg->SetTextAlign(23); latex_seg->SetTextSize(0.045); }
	  latex_seg->DrawLatex(x_center, 1-top-0.015, seg.label);
	}
  }

  if( flag_names && style.description!="" ) {// above the frame
	TLatex *latex_desc = new TLatex();
	latex_desc->SetNDC(); latex_desc->SetTextFont(42); latex_desc->SetTextAlign(11); latex_desc->SetTextSize(0.045);
	latex_desc->DrawLatex(left, 1-top+0.015, style.description);
  }
}

// Plot settings per goodness-of-fit index: the x axes of the target bins (bin index -> variable value, up to two
// segments, e.g. FC and PC), legend positions and labels. Add a case for a new test index; other indices get the
// defaults (bin index axis).
TLeeGoFPlotStyle TLee::Get_GoF_plot_style(int index)
{
  TLeeGoFPlotStyle style;

  switch( index ) {
  case 1:
	style.flag_axis_userAA = 1;
	style.flag_axis_userAB = 1;
	style.title_axis_user = "Reco neutrino energy (MeV)";
	style.axis_user_divisions = 503;
	style.userAA_index_low = 0;  style.userAA_index_hgh = 26; style.userAA_value_low = 0; style.userAA_value_hgh = 2600;
	style.userAB_index_low = 26; style.userAB_index_hgh = 52; style.userAB_value_low = 0; style.userAB_value_hgh = 2600;
	style.legend_x1 = 0.2; style.legend_x2 = 0.4;
	break;

  case 2:
  case 3:
  case 4:
	style.flag_axis_userAA = 1;
	style.title_axis_user = "Reco kinetic energy of #pi^{0} (MeV)";
	style.axis_user_divisions = 508;
	style.userAA_index_low = 0; style.userAA_index_hgh = 11; style.userAA_value_low = 0; style.userAA_value_hgh = 1100;
	break;

  case 5:
  case 8:
	style.flag_axis_userAA = 1;
	style.title_axis_user = "Reco neutrino energy (MeV)";
	style.axis_user_divisions = 508;
	style.userAA_index_low = 0; style.userAA_index_hgh = 26; style.userAA_value_low = 0; style.userAA_value_hgh = 2600;
	break;

  case 6:
	style.flag_axis_userAA = 1;
	style.title_axis_user = "Reco neutrino energy (MeV)";
	style.axis_user_divisions = 508;
	style.userAA_index_low = 0; style.userAA_index_hgh = 18; style.userAA_value_low = 800; style.userAA_value_hgh = 2600;
	break;

  case 7:
	style.flag_axis_userAA = 1;
	style.title_axis_user = "Reco neutrino energy (MeV)";
	style.axis_user_divisions = 508;
	style.userAA_index_low = 0; style.userAA_index_hgh = 8; style.userAA_value_low = 0; style.userAA_value_hgh = 800;
	style.legend_x1 = 0.2; style.legend_x2 = 0.4;
	break;

  case 9:
	style.flag_axis_userAA = 1;
	style.flag_axis_userAB = 1;
	style.title_axis_user = "Reco neutrino energy (MeV)";
	style.axis_user_divisions = 504;
	style.userAA_index_low = 0;  style.userAA_index_hgh = 31; style.userAA_value_low = 0; style.userAA_value_hgh = 3100;
	style.userAB_index_low = 31; style.userAB_index_hgh = 62; style.userAB_value_low = 0; style.userAB_value_hgh = 3100;
	break;

  case 1001:
  case 1002:
	style.legend_total_x1 = 0.2; style.legend_total_x2 = 0.4;
	style.lee_legend = "LEEx";
	style.total_xtitle = "Reco neutrino energy (x100 MeV)";
	break;

  case 1003:
	style.legend_total_x1 = 0.2; style.legend_total_x2 = 0.4;
	style.lee_legend = "Best-fit LEEx";
	style.total_xtitle = "Reco neutrino energy (x100 MeV)";
	break;

  default:
	get_val_plot_style(index, style);// model-validation tests (src/PlotTLee_val.h)
	break;
  }

  return style;
}

// x axis segments of the test index in the current pad (see draw_gof_segments)
void TLee::Draw_GoF_segments(int index, bool flag_axis_labels, double tick_size, double label_size, bool flag_names)
{
  draw_gof_segments(Get_GoF_plot_style(index), flag_axis_labels, tick_size, label_size, flag_names);
}

// Plots of a goodness-of-fit test (Exe_Goodness_of_fit): prediction with uncertainties and data, data/prediction,
// before and, if there are constraining bins (num_X > 0), after the constraint; also file_hists_<index>.root.
void TLee::Plotting_GoF(int index, int num_Y, int num_X, TMatrixD matrix_pred_Y, TMatrixD matrix_data_Y, TMatrixD matrix_YY,
                        double val_chi2_noConstraint, TMatrixD matrix_Y_under_X, TMatrixD matrix_YY_under_XX, double val_chi2_wiConstraint)
{
  TString roostr = "";

  int color_no = kRed;
  int color_wi = kBlue;
  int color_data = kBlack;

  //////////////////////////////////// settings of this test index (Get_GoF_plot_style)

  TLeeGoFPlotStyle style = Get_GoF_plot_style(index);

  bool flag_axis_userAA = style.flag_axis_userAA;
  bool flag_axis_userAB = style.flag_axis_userAB;
  bool flag_axis_user = flag_axis_userAA || flag_axis_userAB || !style.segments.empty();
  int axis_user_divisions = style.axis_user_divisions;
  TString title_axis_user = style.title_axis_user;

  double userAA_index_low = style.userAA_index_low;
  double userAA_index_hgh = style.userAA_index_hgh;
  double userAA_value_low = style.userAA_value_low;
  double userAA_value_hgh = style.userAA_value_hgh;
  double userAA_TickSize       = 0.06;
  double userAA_LabelSize      = 0.078;
  double userAA_clone_TickSize = 0.05;
  double userAA_wi2no_TickSize = 0.03;

  double userAB_index_low = style.userAB_index_low;
  double userAB_index_hgh = style.userAB_index_hgh;
  double userAB_value_low = style.userAB_value_low;
  double userAB_value_hgh = style.userAB_value_hgh;
  double userAB_TickSize       = 0.06;
  double userAB_LabelSize      = 0.078;
  double userAB_clone_TickSize = 0.05;
  double userAB_wi2no_TickSize = 0.03;

  TLine *line_FC_PC = new TLine(num_Y/2, 0, num_Y/2, 2);
  line_FC_PC->SetLineColor(kGray+1);
  line_FC_PC->SetLineWidth(4);

  ///////////

  TGaxis *axis_userAA = new TGaxis(userAA_index_low, 0, userAA_index_hgh, 0,   userAA_value_low, userAA_value_hgh, axis_user_divisions, "S");
  axis_userAA->SetName("axis_userAA");
  axis_userAA->SetTickSize(userAA_TickSize);
  axis_userAA->SetLabelSize(userAA_LabelSize);
  axis_userAA->SetLabelFont(42);
  TGaxis *axis_userAA_clone = (TGaxis*)axis_userAA->Clone("axis_userAA_clone");
  axis_userAA_clone->SetLabelSize(0);
  axis_userAA_clone->SetTickSize(userAA_clone_TickSize);
  TGaxis *axis_userAA_wi2no = (TGaxis*)axis_userAA->Clone("axis_userAA_clone");
  axis_userAA_wi2no->SetLabelSize(userAA_clone_TickSize);
  axis_userAA_wi2no->SetTickSize(userAA_wi2no_TickSize);

  TGaxis *axis_userAB = new TGaxis(userAB_index_low, 0, userAB_index_hgh, 0,   userAB_value_low, userAB_value_hgh, axis_user_divisions, "S");
  axis_userAB->SetName("axis_userAB");
  axis_userAB->SetTickSize(userAB_TickSize);
  axis_userAB->SetLabelSize(userAB_LabelSize);
  axis_userAB->SetLabelFont(42);
  TGaxis *axis_userAB_clone = (TGaxis*)axis_userAB->Clone("axis_userAB_clone");
  axis_userAB_clone->SetLabelSize(0);
  axis_userAB_clone->SetTickSize(userAB_clone_TickSize);
  TGaxis *axis_userAB_wi2no = (TGaxis*)axis_userAB->Clone("axis_userAB_clone");
  axis_userAB_wi2no->SetLabelSize(userAB_clone_TickSize);
  axis_userAB_wi2no->SetTickSize(userAB_wi2no_TickSize);


  ///////////////////////////////////////////////////////////////////////////////////////////// noConstraint

  roostr = TString::Format("h1_pred_Y_noConstraint_%02d", index);
  TH1D *h1_pred_Y_noConstraint = new TH1D(roostr, "", num_Y, 0, num_Y);
  for(int ibin=1; ibin<=num_Y; ibin++) {
	h1_pred_Y_noConstraint->SetBinContent( ibin, matrix_pred_Y(ibin-1, 0) );
	double val_err = sqrt( matrix_YY(ibin-1, ibin-1) );
	h1_pred_Y_noConstraint->SetBinError( ibin, val_err );
  }

  TGraphAsymmErrors *gh_data = new TGraphAsymmErrors();
  TGraphAsymmErrors *gh_ratio_noConstraint = new TGraphAsymmErrors();
  map<int, double>array_val_data_low;
  map<int, double>array_val_data_hgh;

  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_data = matrix_data_Y(ibin-1, 0);
	double val_pred_noConstraint = matrix_pred_Y(ibin-1, 0);

	double val_data_low = 0;
	double val_data_hgh = 0;
	int idx_data = (int)(val_data+0.5);
	if( idx_data>100 ) {
	  val_data_low = val_data - sqrt(val_data);
	  val_data_hgh = val_data + sqrt(val_data);
	}
	else {
	  val_data_low = DataBase::yl[ idx_data ];
	  val_data_hgh = DataBase::yh[ idx_data ];
	}
	gh_data->SetPoint( ibin-1, ibin-0.5, val_data );
	gh_data->SetPointError( ibin-1, 0.5, 0.5, val_data-val_data_low, val_data_hgh-val_data );

	array_val_data_low[ibin-1] = val_data_low;
	array_val_data_hgh[ibin-1] = val_data_hgh;

	double val_ratio_no = val_data/val_pred_noConstraint;
	double val_ratio_no_low = val_ratio_no - val_data_low/val_pred_noConstraint;
	double val_ratio_no_hgh = val_data_hgh/val_pred_noConstraint - val_ratio_no;
	if( val_ratio_no!=val_ratio_no || std::isinf(val_ratio_no) ) val_ratio_no = 0;
	gh_ratio_noConstraint->SetPoint( ibin-1, ibin-0.5, val_ratio_no );
	gh_ratio_noConstraint->SetPointError( ibin-1, 0.5, 0.5, val_ratio_no_low, val_ratio_no_hgh );
  }

  ///////

  double ymax_pred = 0;
  double ymax_data = 0;

  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_pred = h1_pred_Y_noConstraint->GetBinContent(ibin)
	                  + h1_pred_Y_noConstraint->GetBinError(ibin);
	if( ymax_pred<val_pred ) ymax_pred = val_pred;

	double xx_data(0), yy_data(0);
	gh_data->GetPoint(ibin-1, xx_data, yy_data);
	if( ymax_data<yy_data ) ymax_data = yy_data;
  }

  ///////

  roostr = TString::Format("canv_spectra_GoF_no_%02d", index);
  TCanvas *canv_spectra_GoF_no = new TCanvas(roostr, roostr, 1000, 950);

  ///////
  canv_spectra_GoF_no->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.05);
  pad_top_no->Draw(); pad_top_no->cd();

  TH1D *h1_pred_Y_noConstraint_clone = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_clone");

  h1_pred_Y_noConstraint->Draw("e2");
  h1_pred_Y_noConstraint->SetMinimum(0.);
  //h1_pred_Y_noConstraint->SetMaximum(4.6);
  if( ymax_data>ymax_pred/**1.05*/ ) h1_pred_Y_noConstraint->SetMaximum(ymax_data*1.1);
  h1_pred_Y_noConstraint->SetMarkerSize(0.);
  h1_pred_Y_noConstraint->SetFillColor(color_no); h1_pred_Y_noConstraint->SetFillStyle(3005);
  h1_pred_Y_noConstraint->SetLineColor(color_no);
  func_title_size(h1_pred_Y_noConstraint, 0.065, 0.065, 0.065, 0.065);
  func_xy_title(h1_pred_Y_noConstraint, "", "Entries");
  h1_pred_Y_noConstraint->GetXaxis()->SetLabelOffset(2);
  h1_pred_Y_noConstraint->GetYaxis()->CenterTitle();
  h1_pred_Y_noConstraint->GetYaxis()->SetTitleOffset(1.2);
  h1_pred_Y_noConstraint->GetYaxis()->SetTickLength(0.02);

  h1_pred_Y_noConstraint_clone->Draw("same hist");
  h1_pred_Y_noConstraint_clone->SetLineColor(color_no);

  gh_data->Draw("same pe");
  gh_data->SetMarkerStyle(20); gh_data->SetMarkerSize(1.12);
  gh_data->SetMarkerColor(color_data); gh_data->SetLineColor(color_data);
  if( num_X==0 ) {
	gh_data->SetMarkerColor(kBlue); gh_data->SetLineColor(kBlue);
  }

  h1_pred_Y_noConstraint->Draw("same axis");

  TLegend *lg_top_no = new TLegend(0.5, 0.60, 0.85, 0.85);
  if (moveleg){
    //lg_top_total->SetX1(0.25); lg_top_total->SetX2(0.6);
    lg_top_no->SetX1(0.2); lg_top_no->SetX2(0.55);
  }
  if( style.legend_x1>=0 ) { lg_top_no->SetX1(style.legend_x1); lg_top_no->SetX2(style.legend_x2);}
  lg_top_no->AddEntry(gh_data, "Data", "lep");
  lg_top_no->AddEntry(h1_pred_Y_noConstraint, TString::Format("#color[%d]{Pred no constraint}", color_no), "lf");
  lg_top_no->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_no, val_chi2_noConstraint, num_Y), "");
  lg_top_no->Draw();
  lg_top_no->SetBorderSize(0); lg_top_no->SetFillStyle(0); lg_top_no->SetTextSize(0.065);

  ///////
  canv_spectra_GoF_no->cd();
  TPad *pad_bot_no = new TPad("pad_bot_no", "pad_bot_no", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_no->Draw(); pad_bot_no->cd();

  TH1D *h1_pred_Y_noConstraint_rel_error = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_rel_error");
  h1_pred_Y_noConstraint_rel_error->Reset();
  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin);
	double val_err = h1_pred_Y_noConstraint->GetBinError(ibin);
	double rel_err = val_err/val_cv;
	if( val_cv==0 ) rel_err = 0;
	h1_pred_Y_noConstraint_rel_error->SetBinContent(ibin, 1);
	h1_pred_Y_noConstraint_rel_error->SetBinError(ibin, rel_err);
  }

  h1_pred_Y_noConstraint_rel_error->Draw("e2");
  h1_pred_Y_noConstraint_rel_error->SetMinimum(0); h1_pred_Y_noConstraint_rel_error->SetMaximum(2);
  func_title_size(h1_pred_Y_noConstraint_rel_error, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_pred_Y_noConstraint_rel_error, "Bin index", "Data / Pred");
  h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetTickLength(0.05);
  h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetLabelOffset(0.005);
  h1_pred_Y_noConstraint_rel_error->GetXaxis()->CenterTitle(); h1_pred_Y_noConstraint_rel_error->GetYaxis()->CenterTitle();
  h1_pred_Y_noConstraint_rel_error->GetYaxis()->SetTitleOffset(0.99);
  h1_pred_Y_noConstraint_rel_error->GetYaxis()->SetNdivisions(509);

  TF1 *line_no = new TF1("line_no", "1", 0, 1e6); line_no->Draw("same");
  line_no->SetLineColor(kBlack); line_no->SetLineStyle(7);

  gh_ratio_noConstraint->Draw("same pe");
  gh_ratio_noConstraint->SetMarkerStyle(20); gh_ratio_noConstraint->SetMarkerSize(1.12);
  gh_ratio_noConstraint->SetMarkerColor(color_no); gh_ratio_noConstraint->SetLineColor(color_no);

  h1_pred_Y_noConstraint_rel_error->Draw("same axis");

  if( flag_axis_user ) {
	///////////////////// bot
	h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetTickLength(0);
	h1_pred_Y_noConstraint_rel_error->GetXaxis()->SetLabelSize(0);
	h1_pred_Y_noConstraint_rel_error->SetXTitle( title_axis_user );

	if( flag_axis_userAA && flag_axis_userAB ) line_FC_PC->Draw("same");

	if( flag_axis_userAA ) axis_userAA->Draw();
	if( flag_axis_userAB ) axis_userAB->Draw();
	draw_gof_segments(style, true, 0.06, 0.06, false);

	///////////////////// top
	canv_spectra_GoF_no->cd(); pad_top_no->cd();
	h1_pred_Y_noConstraint->GetXaxis()->SetTickLength(0);
	h1_pred_Y_noConstraint->GetXaxis()->SetLabelSize(0);
	if( flag_axis_userAA ) axis_userAA_clone->Draw();
	if( flag_axis_userAB ) axis_userAB_clone->Draw();
	draw_gof_segments(style, false, 0.05, 0, true);
  }

  if( num_X==0 ) {
	gh_ratio_noConstraint->SetMarkerColor(kBlue); gh_ratio_noConstraint->SetLineColor(kBlue);

	roostr = TString::Format("canv_spectra_GoF_no_%02d.png", index);
	canv_spectra_GoF_no->SaveAs(roostr);
	return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////// wiConstraint

  roostr = TString::Format("h1_pred_Y_wiConstraint_%02d", index);
  TH1D *h1_pred_Y_wiConstraint = (TH1D*)h1_pred_Y_noConstraint->Clone(roostr);
  for(int ibin=1; ibin<=num_Y; ibin++) {
	h1_pred_Y_wiConstraint->SetBinContent( ibin, matrix_Y_under_X(ibin-1, 0) );
	double val_err = sqrt( matrix_YY_under_XX(ibin-1, ibin-1) );
	h1_pred_Y_wiConstraint->SetBinError( ibin, val_err );
  }

  TGraphAsymmErrors *gh_ratio_wiConstraint = new TGraphAsymmErrors();

  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_data = matrix_data_Y(ibin-1, 0);
	double val_pred_wiConstraint = matrix_Y_under_X(ibin-1, 0);

	double val_data_low = array_val_data_low[ibin-1];
	double val_data_hgh = array_val_data_hgh[ibin-1];
	gh_data->SetPoint( ibin-1, ibin-0.5, val_data );
	gh_data->SetPointError( ibin-1, 0.5, 0.5, val_data-val_data_low, val_data_hgh-val_data );

	///////

	double val_ratio_wi = val_data/val_pred_wiConstraint;
	double val_ratio_wi_low = val_ratio_wi - val_data_low/val_pred_wiConstraint;
	double val_ratio_wi_hgh = val_data_hgh/val_pred_wiConstraint - val_ratio_wi;
	if( val_ratio_wi!=val_ratio_wi || std::isinf(val_ratio_wi) ) val_ratio_wi = 0;
	gh_ratio_wiConstraint->SetPoint( ibin-1, ibin-0.5, val_ratio_wi );
	gh_ratio_wiConstraint->SetPointError( ibin-1, 0.5, 0.5, val_ratio_wi_low, val_ratio_wi_hgh );
  }

  roostr = TString::Format("canv_spectra_GoF_wi_%02d", index);
  TCanvas *canv_spectra_GoF_wi = new TCanvas(roostr, roostr, 1000, 950);

  ///////
  canv_spectra_GoF_wi->cd();
  TPad *pad_top_wi = new TPad("pad_top_wi", "pad_top_wi", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_wi, 0.15, 0.1, 0.1, 0.05);
  pad_top_wi->Draw(); pad_top_wi->cd();

  TH1D *h1_pred_Y_wiConstraint_clone = (TH1D*)h1_pred_Y_wiConstraint->Clone("h1_pred_Y_wiConstraint_clone");

  h1_pred_Y_wiConstraint->Draw("e2");
  h1_pred_Y_wiConstraint->SetFillColor(color_wi); h1_pred_Y_wiConstraint->SetFillStyle(3004);
  h1_pred_Y_wiConstraint->SetLineColor(color_wi);

  h1_pred_Y_wiConstraint_clone->Draw("same hist");
  h1_pred_Y_wiConstraint_clone->SetLineColor(color_wi); h1_pred_Y_wiConstraint_clone->SetFillStyle(0);

  gh_data->Draw("same pe");

  h1_pred_Y_wiConstraint->Draw("same axis");

  TLegend *lg_top_wi = new TLegend(0.5, 0.57, 0.85, 0.87);
  if (moveleg){
    //lg_top_total->SetX1(0.25); lg_top_total->SetX2(0.6);
    lg_top_wi->SetX1(0.2); lg_top_wi->SetX2(0.55);
  }
  if( style.legend_x1>=0 ) { lg_top_wi->SetX1(style.legend_x1); lg_top_wi->SetX2(style.legend_x2);}
  lg_top_wi->AddEntry(gh_data, "Data", "lep");
  lg_top_wi->AddEntry(h1_pred_Y_wiConstraint, TString::Format("#color[%d]{Pred with constraint}", color_wi), "lf");
  lg_top_wi->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_wi, val_chi2_wiConstraint, num_Y), "");
  lg_top_wi->AddEntry("", TString::Format("#color[%d]{p-value: %3.2f}", color_wi, TMath::Prob(val_chi2_wiConstraint, num_Y)), "");
  lg_top_wi->AddEntry("", TString::Format("#color[%d]{%3.2f#sigma}", color_wi, RooStats::PValueToSignificance(TMath::Prob(val_chi2_wiConstraint, num_Y) / 2.0)), "");
  lg_top_wi->Draw();
  lg_top_wi->SetBorderSize(0); lg_top_wi->SetFillStyle(0); lg_top_wi->SetTextSize(0.065);

  ///////
  canv_spectra_GoF_wi->cd();
  TPad *pad_bot_wi = new TPad("pad_bot_wi", "pad_bot_wi", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_wi, 0.15, 0.1, 0.05, 0.3);
  pad_bot_wi->Draw(); pad_bot_wi->cd();

  TH1D *h1_pred_Y_wiConstraint_rel_error = (TH1D*)h1_pred_Y_noConstraint_rel_error->Clone("h1_pred_Y_wiConstraint_rel_error");
  h1_pred_Y_wiConstraint_rel_error->Reset();
  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin);
	double val_err = h1_pred_Y_wiConstraint->GetBinError(ibin);
	double rel_err = val_err/val_cv;
	if( val_cv==0 ) rel_err = 0;
	h1_pred_Y_wiConstraint_rel_error->SetBinContent(ibin, 1);
	h1_pred_Y_wiConstraint_rel_error->SetBinError(ibin, rel_err);

  }

  h1_pred_Y_wiConstraint_rel_error->Draw("e2");
  h1_pred_Y_wiConstraint_rel_error->SetFillColor(color_wi); h1_pred_Y_wiConstraint_rel_error->SetFillStyle(3004);

  TF1 *line_wi = new TF1("line_wi", "1", 0, 1e6); line_wi->Draw("same");
  line_wi->SetLineColor(kBlack); line_wi->SetLineStyle(7);

  gh_ratio_wiConstraint->Draw("same pe");
  gh_ratio_wiConstraint->SetMarkerStyle(20); gh_ratio_wiConstraint->SetMarkerSize(1.12);
  gh_ratio_wiConstraint->SetMarkerColor(color_wi); gh_ratio_wiConstraint->SetLineColor(color_wi);

  h1_pred_Y_wiConstraint_rel_error->Draw("same axis");

  if( flag_axis_user ) {
	///////////////////// bot

	if( flag_axis_userAA && flag_axis_userAB ) line_FC_PC->Draw("same");

	if( flag_axis_userAA ) axis_userAA->Draw();
	if( flag_axis_userAB ) axis_userAB->Draw();
	draw_gof_segments(style, true, 0.06, 0.06, false);

	///////////////////// top
	canv_spectra_GoF_wi->cd(); pad_top_wi->cd();
	if( flag_axis_userAA ) axis_userAA_clone->Draw();
	if( flag_axis_userAB ) axis_userAB_clone->Draw();
	draw_gof_segments(style, false, 0.05, 0, true);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////

  canv_spectra_GoF_no->cd(); pad_top_no->cd(); pad_top_no->Update(); double y2_no = gPad->GetUymax();
  canv_spectra_GoF_wi->cd(); pad_top_wi->cd(); pad_top_wi->Update(); double y2_wi = gPad->GetUymax();

  if( y2_no > y2_wi ) {
	canv_spectra_GoF_wi->cd(); pad_top_wi->cd(); h1_pred_Y_wiConstraint->SetMaximum(y2_no);
	h1_pred_Y_wiConstraint->Draw("same e2"); h1_pred_Y_wiConstraint_clone->Draw("same hist");
	gh_data->Draw("same pe"); h1_pred_Y_wiConstraint->Draw("same axis");
  }
  else {
	canv_spectra_GoF_no->cd(); pad_top_no->cd(); h1_pred_Y_noConstraint->SetMaximum(y2_wi);
	h1_pred_Y_noConstraint->Draw("same e2"); h1_pred_Y_noConstraint_clone->Draw("same hist");
	gh_data->Draw("same pe"); h1_pred_Y_noConstraint->Draw("same axis");
  }

  roostr = TString::Format("canv_spectra_GoF_no_%02d.png", index); canv_spectra_GoF_no->SaveAs(roostr);
  roostr = TString::Format("canv_spectra_GoF_wi_%02d.png", index); canv_spectra_GoF_wi->SaveAs(roostr);

  ////////////////

  roostr = TString::Format("h1_spectra_wi2no_%02d", index);
  TString roostr_wi2no = roostr;
  //TH1D *h1_spectra_wi2no = (TH1D*)h1_pred_Y_noConstraint->Clone(roostr);
  TH1D *h1_spectra_wi2no = new TH1D(roostr, "", num_Y, 0, num_Y);
  h1_spectra_wi2no->SetFillStyle(0);
  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_wiConstraint = h1_pred_Y_wiConstraint->GetBinContent(ibin);
	double val_noConstraint = h1_pred_Y_noConstraint->GetBinContent(ibin);
	double val_wi2no = val_wiConstraint/val_noConstraint;
	if( val_noConstraint==0 ) val_wi2no = 0;
	h1_spectra_wi2no->SetBinContent(ibin, val_wi2no);
  }

  roostr = TString::Format("canv_spectra_wi2no_%02d", index);
  TCanvas *canv_spectra_wi2no = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_spectra_wi2no, 0.15, 0.1, 0.1, 0.15);
  canv_spectra_wi2no->SetGridy();

  h1_spectra_wi2no->Draw("hist");
  h1_spectra_wi2no->SetLineColor(color_wi);

  h1_spectra_wi2no->SetMinimum(0);
  h1_spectra_wi2no->SetMaximum(2);
  func_title_size(h1_spectra_wi2no, 0.05, 0.05, 0.05, 0.05);

  func_xy_title(h1_spectra_wi2no, "Bin index", "Prediction wi/no constraint");
  h1_spectra_wi2no->GetXaxis()->CenterTitle(); h1_spectra_wi2no->GetYaxis()->CenterTitle();
  h1_spectra_wi2no->GetYaxis()->SetTitleOffset(1.18);
  h1_spectra_wi2no->GetYaxis()->SetNdivisions(509);
  h1_spectra_wi2no->GetYaxis()->SetTickLength(0.03);
  if( flag_axis_user ) {/// ttt
	func_xy_title(h1_spectra_wi2no, title_axis_user, "Prediction wi/no constraint");
	h1_spectra_wi2no->GetXaxis()->SetLabelSize(0);
	if( flag_axis_userAA ) axis_userAA_wi2no->Draw();
	if( flag_axis_userAB ) axis_userAB_wi2no->Draw();
	draw_gof_segments(style, true, 0.03, 0.04, true);
  }

  //roostr = TString::Format("canv_spectra_wi2no_%02d.png", index); canv_spectra_wi2no->SaveAs(roostr);

  //h1_spectra_wi2no->SaveAs("file_h1_spectra_wi2no.root");

  // for(int ibin=1; ibin<=8; ibin++) {
  //   double cv_no = h1_pred_Y_noConstraint->GetBinContent(ibin);
  //   double err_no = h1_pred_Y_noConstraint->GetBinError(ibin);
  //   double cv_wi = h1_pred_Y_wiConstraint->GetBinContent(ibin);
  //   double err_wi = h1_pred_Y_wiConstraint->GetBinError(ibin);

  //   cout<<TString::Format(" ---> %d, (no con) %5.2f %5.2f, relerr %5.2f, (wi con) %5.2f %5.2f, relerr %5.2f",
  // 			  ibin, cv_no, err_no, err_no/cv_no, cv_wi, err_wi, err_wi/cv_wi)<<endl;
  // }

  ////////////////

  roostr = TString::Format("canv_spectra_GoF_total_%02d", index);
  TCanvas *canv_spectra_GoF_total = new TCanvas(roostr, roostr, 1000, 950);

  ///////
  canv_spectra_GoF_total->cd();
  TPad *pad_top_total = new TPad("pad_top_total", "pad_top_total", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_total, 0.15, 0.1, 0.1, 0.05);
  pad_top_total->Draw(); pad_top_total->cd();

  h1_pred_Y_wiConstraint->Draw("e2");
  h1_pred_Y_wiConstraint->SetXTitle("");// top pad: the x axis is labelled in the bottom pad
  //if( index==7 ) h1_pred_Y_wiConstraint->SetMaximum(25);
  //if( index==9 ) h1_pred_Y_wiConstraint->SetMaximum(50);
  h1_pred_Y_noConstraint->Draw("same e2");
  h1_pred_Y_wiConstraint_clone->Draw("same hist");
  h1_pred_Y_noConstraint_clone->Draw("same hist");
  gh_data->Draw("same pe");
  h1_pred_Y_wiConstraint->Draw("same axis");

  /*
     cout<<endl;
     double data_FC = 0;
     double data_PC = 0;
     for( int i=0; i<8; i++ ) {
     data_FC += matrix_data_Y(i, 0);
     data_PC += matrix_data_Y(i+8, 0);
     }
     cout<<" ---> data "<< data_FC<<"\t"<<data_PC<<endl;
     cout<<" ---> pred noConstraint "<<h1_pred_Y_noConstraint->Integral(1,8)<<"\t"<<h1_pred_Y_noConstraint->Integral(9,16)<<endl;
     cout<<" ---> pred wiConstraint "<<h1_pred_Y_wiConstraint->Integral(1,8)<<"\t"<<h1_pred_Y_wiConstraint->Integral(9,16)<<endl;
     cout<<endl;
   */

  TLegend *lg_top_total = new TLegend(0.5, 0.45, 0.85, 0.85);
  if (moveleg){
    //lg_top_total->SetX1(0.25); lg_top_total->SetX2(0.6);
    lg_top_total->SetX1(0.2); lg_top_total->SetX2(0.55);
  }
  // h1_pred_Y_wiConstraint->SetMaximum(40);
  // lg_top_total->SetX1(0.55); lg_top_total->SetX2(0.95);
  // lg_top_total->SetX1(0.2); lg_top_total->SetX2(0.4);
  if( style.legend_x1>=0 ) { lg_top_total->SetX1(style.legend_x1); lg_top_total->SetX2(style.legend_x2);}
  if( style.legend_total_x1>=0 ) { lg_top_total->SetX1(style.legend_total_x1); lg_top_total->SetX2(style.legend_total_x2);}
  if( style.lee_legend!="" ) lg_top_total->AddEntry("", TString::Format("#color[%d]{%s = %3.1f}", kGreen+1, style.lee_legend.Data(), scaleF_Lee), "");

  lg_top_total->AddEntry(gh_data, "Data", "lep");
  lg_top_total->AddEntry(h1_pred_Y_noConstraint, TString::Format("#color[%d]{Pred no constraint}", color_no), "lf");
  lg_top_total->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_no, val_chi2_noConstraint, num_Y), "");
  lg_top_total->AddEntry(h1_pred_Y_wiConstraint, TString::Format("#color[%d]{Pred with constraint}", color_wi), "lf");
  lg_top_total->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_wi, val_chi2_wiConstraint, num_Y), "");
  lg_top_total->Draw();
  lg_top_total->SetBorderSize(0); lg_top_total->SetFillStyle(0); lg_top_total->SetTextSize(0.065);

  ///////
  canv_spectra_GoF_total->cd();
  TPad *pad_bot_total = new TPad("pad_bot_total", "pad_bot_total", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_total, 0.15, 0.1, 0.05, 0.3);
  pad_bot_total->Draw(); pad_bot_total->cd();
  //pad_bot_total->SetTicky();

  h1_pred_Y_noConstraint_rel_error->Draw("e2");
  //h1_pred_Y_noConstraint_rel_error->SetYTitle("Uncertainties");
  h1_pred_Y_wiConstraint_rel_error->Draw("same e2");

  TF1 *line_total = new TF1("line_total", "1", 0, 1e6); line_total->Draw("same");
  line_total->SetLineColor(kBlack); line_total->SetLineStyle(7);

  gh_ratio_noConstraint->Draw("same pe");
  gh_ratio_wiConstraint->Draw("same pe");

  if( flag_axis_user ) {
	///////////////////// bot

	if( flag_axis_userAA && flag_axis_userAB ) line_FC_PC->Draw("same");

	if( flag_axis_userAA ) axis_userAA->Draw();
	if( flag_axis_userAB ) axis_userAB->Draw();
	draw_gof_segments(style, true, 0.06, 0.06, false);

	///////////////////// top
	canv_spectra_GoF_total->cd(); pad_top_total->cd();
	if( flag_axis_userAA ) axis_userAA_clone->Draw();
	if( flag_axis_userAB ) axis_userAB_clone->Draw();
	draw_gof_segments(style, false, 0.05, 0, true);
  }

  if( style.total_xtitle!="" ) h1_pred_Y_noConstraint_rel_error->SetXTitle( style.total_xtitle );

  // h1_spectra_wi2no->Draw("same");
  // h1_spectra_wi2no->SetLineColor(kGreen+1);
  // TLegend *lg_wi2no = new TLegend(0.92, 0.15, 0.94, 0.60);
  // lg_wi2no->SetHeader( TString::Format("#color[%d]{Prediction wi/wo}", kGreen+1) );
  // lg_wi2no->Draw("same"); lg_wi2no->SetTextSize(0.078); lg_wi2no->SetTextAngle(90);
  // lg_wi2no->SetBorderSize(0);

  // TLatex *latex = new TLatex(0.5, 0.5, TString::Format("#color[%d]{Predictioin wi/wo}", kGreen+1));
  // latex->Draw("same"); latex->SetTextSize(0.078); //latex->SetTextAngle(90);

  // h1_pred_Y_noConstraint_rel_error->Draw("same axis");
  // TLegend *lg_bot_total = new TLegend(0.5, 0.85, 0.85, 0.93);
  // lg_bot_total->AddEntry(h1_pred_Y_noConstraint_rel_error, TString::Format("#color[%d]{Prediction wi/no}", kGreen+1), "l");
  // lg_bot_total->Draw();
  // lg_bot_total->SetBorderSize(0); lg_bot_total->SetTextSize(0.078);
  // lg_bot_total->SetFillColor(10);

  roostr = TString::Format("canv_spectra_GoF_total_%02d.png", index); canv_spectra_GoF_total->SaveAs(roostr);

  //Erin
  TFile *file_hists = new TFile(TString::Format("file_hists_%02d.root", index), "recreate");
  file_hists->cd();
  gh_data->SetName("gh_data");
  gh_data->Write();
  h1_pred_Y_wiConstraint->Write();
  h1_pred_Y_noConstraint->Write();
  h1_pred_Y_wiConstraint_clone->Write();
  h1_pred_Y_noConstraint_clone->Write();
  h1_pred_Y_noConstraint_rel_error->Write();
  h1_pred_Y_wiConstraint_rel_error->Write();
  gh_ratio_noConstraint->SetName("gh_ratio_noConstraint");
  gh_ratio_noConstraint->Write();
  gh_ratio_wiConstraint->SetName("gh_ratio_wiConstraint");
  gh_ratio_wiConstraint->Write();
  matrix_data_Y.Write("matrix_data_Y");
  matrix_Y_under_X.Write("matrix_Y_under_X");
  matrix_YY_under_XX.Write("matrix_YY_under_XX");
  file_hists->Write();
  file_hists->Close();


  //////////////////////////////////////////////////////////////////

  roostr = TString::Format("h1_spectra_relerr_%02d", index);
  TH1D *h1_spectra_relerr = new TH1D(roostr, "", num_Y, 0, num_Y);

  roostr = TString::Format("h1_spectra_relerr_wi_%02d", index);
  TH1D *h1_spectra_relerr_wi = new TH1D(roostr, "", num_Y, 0, num_Y);

  for(int ibin=1; ibin<=num_Y; ibin++) {
	double val_noConstraint = h1_pred_Y_noConstraint_rel_error->GetBinError(ibin);
	double val_wiConstraint = h1_pred_Y_wiConstraint_rel_error->GetBinError(ibin);
	h1_spectra_relerr->SetBinContent(ibin, val_noConstraint);
	h1_spectra_relerr_wi->SetBinContent(ibin, val_wiConstraint);
  }

  roostr = TString::Format("canv_spectra_relerr_%02d", index);
  TCanvas *canv_spectra_relerr = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_spectra_relerr, 0.15, 0.1, 0.1, 0.15);
  canv_spectra_relerr->SetGridy();

  h1_spectra_relerr->Draw("hist");
  h1_spectra_relerr->SetLineColor(color_no);

  h1_spectra_relerr_wi->Draw("same hist");
  h1_spectra_relerr_wi->SetLineColor(color_wi);

  h1_spectra_relerr->SetMinimum(0);
  //h1_spectra_relerr->SetMaximum(0.5);
  func_title_size(h1_spectra_relerr, 0.05, 0.05, 0.05, 0.05);

  h1_spectra_relerr->Draw("same axis");

  func_xy_title(h1_spectra_relerr, "Bin index", "Rel.Err to Pred no constraint");
  h1_spectra_relerr->GetXaxis()->CenterTitle(); h1_spectra_relerr->GetYaxis()->CenterTitle();
  h1_spectra_relerr->GetYaxis()->SetTitleOffset(1.18);
  h1_spectra_relerr->GetYaxis()->SetNdivisions(509);
  h1_spectra_relerr->GetYaxis()->SetTickLength(0.03);
  if( flag_axis_user ) {/// ttt
	func_xy_title(h1_spectra_relerr, title_axis_user,"Rel.Err to Pred no constraint");

	if( flag_axis_userAA && flag_axis_userAB ) {
	  line_FC_PC->Draw("same");
	  line_FC_PC->SetY2(0.5);
	}

	h1_spectra_relerr->GetXaxis()->SetLabelSize(0);
	h1_spectra_relerr->GetXaxis()->SetTickSize(0);

	if( flag_axis_userAA )  {
	  axis_userAA->Draw();
	  axis_userAA->SetTickSize(0.06);
	  axis_userAA->SetLabelSize(0.05);
	}
	if( flag_axis_userAB )  {
	  axis_userAB->Draw();
	  axis_userAB->SetTickSize(0.06);
	  axis_userAB->SetLabelSize(0.05);
	}
	draw_gof_segments(style, true, 0.03, 0.04, true);
  }

  roostr = TString::Format("canv_spectra_relerr_%02d.png", index); canv_spectra_relerr->SaveAs(roostr);
  //roostr = TString::Format("canv_h1_spectra_relerr_wi_%02d.root", index); h1_spectra_relerr_wi->SaveAs(roostr);

}
