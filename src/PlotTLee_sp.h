// Plot settings of the muon spacepoint data/MC comparisons (numuCCana/configurations/opendata2/sp*_val/TLee_config.txt),
// by test index: the x axis of the target bins (one segment per target channel) and a description of the test. Generated
// by numucc/studies/make_sp_configs.py, do not edit by hand. Included by src/PlotTLee.h (Get_GoF_plot_style).

static bool get_sp_plot_style(int index, TLeeGoFPlotStyle& style)
{
  switch( index ) {
  case 60001:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 60002:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 60003:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60004:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60005:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60006:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60007:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60008:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60009:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60010:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60011:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 60012:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 60013:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 60014:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 60015:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60016:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60017:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60018:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60019:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60020:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 60021:   // sp11_inc_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} LFC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 60022:   // sp11_inc_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} LPC (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 60023:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, "LFC"}, {31, 61, 0, 150000, "LPC"} };
	return true;
  case 60024:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60025:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60026:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60027:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60028:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, "LFC"}, {31, 61, 0, 150000, "LPC"} };
	return true;
  case 60029:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, "LFC"}, {31, 61, 0, 120000, "LPC"} };
	return true;
  case 60030:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60031:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60032:   // sp11_inc_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "LFC"}, {26, 51, 0, 100000, "LPC"} };
	return true;
  case 60033:   // sp11_inc_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med}, all channels (sp11_inc_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, "LFC"}, {29, 57, 18000, 46000, "LPC"} };
	return true;
  case 61001:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61002:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61003:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61004:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61005:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61006:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61007:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61008:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61009:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61010:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61011:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61012:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61013:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61014:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61015:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61016:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61017:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61018:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61019:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61020:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61021:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61022:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61023:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61024:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 61025:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 61026:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 61027:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 61028:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 61029:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61030:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61031:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61032:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61033:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61034:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61035:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61036:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61037:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61038:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61039:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61040:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 61041:   // sp11_wc0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} Np LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 61042:   // sp11_wc0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} Np LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 61043:   // sp11_wc0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} WC0p LFC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 61044:   // sp11_wc0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} WC0p LPC (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 61045:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, "Np LFC"}, {31, 61, 0, 150000, "Np LPC"}, {62, 92, 0, 150000, "WC0p LFC"}, {93, 123, 0, 150000, "WC0p LPC"} };
	return true;
  case 61046:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61047:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61048:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61049:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61050:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, "Np LFC"}, {31, 61, 0, 150000, "Np LPC"}, {62, 92, 0, 150000, "WC0p LFC"}, {93, 123, 0, 150000, "WC0p LPC"} };
	return true;
  case 61051:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, "Np LFC"}, {31, 61, 0, 120000, "Np LPC"}, {62, 92, 0, 120000, "WC0p LFC"}, {93, 123, 0, 120000, "WC0p LPC"} };
	return true;
  case 61052:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61053:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61054:   // sp11_wc0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "WC0p LFC"}, {78, 103, 0, 100000, "WC0p LPC"} };
	return true;
  case 61055:   // sp11_wc0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med}, all channels (sp11_wc0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, "Np LFC"}, {29, 57, 18000, 46000, "Np LPC"}, {58, 86, 18000, 46000, "WC0p LFC"}, {87, 115, 18000, 46000, "WC0p LPC"} };
	return true;
  case 62001:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62002:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62003:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62004:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62005:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62006:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62007:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62008:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62009:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62010:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62011:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62012:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62013:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62014:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62015:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62016:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62017:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62018:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62019:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62020:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62021:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62022:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62023:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62024:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62025:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62026:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62027:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62028:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62029:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62030:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62031:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62032:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62033:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62034:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62035:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62036:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62037:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62038:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62039:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62040:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62041:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62042:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62043:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62044:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62045:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62046:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62047:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62048:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 62049:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62050:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62051:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62052:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62053:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62054:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62055:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62056:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 120000, ""} };
	return true;
  case 62057:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62058:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62059:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62060:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62061:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62062:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62063:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62064:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62065:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62066:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62067:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62068:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62069:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62070:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62071:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62072:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62073:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62074:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62075:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62076:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62077:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62078:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62079:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62080:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 25, 0, 100000, ""} };
	return true;
  case 62081:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} Np LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62082:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} Np LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62083:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} PorLNp LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62084:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} PorLNp LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62085:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} WCTrk0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62086:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} WCTrk0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62087:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} WCG0p LFC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62088:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med} WCG0p LPC (sp11_wcx0pNp_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 18000, 46000, ""} };
	return true;
  case 62089:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the vertex [(q+1000)#times10]";
	style.description = "q_{0}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 30, 0, 150000, "Np LFC"}, {31, 61, 0, 150000, "Np LPC"}, {62, 92, 0, 150000, "PorLNp LFC"}, {93, 123, 0, 150000, "PorLNp LPC"}, {124, 154, 0, 150000, "WCTrk0p LFC"}, {155, 185, 0, 150000, "WCTrk0p LPC"}, {186, 216, 0, 150000, "WCG0p LFC"}, {217, 247, 0, 150000, "WCG0p LPC"} };
	return true;
  case 62090:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the vertex [(q+1000)#times10]";
	style.description = "q_{1}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62091:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the vertex [(q+1000)#times10]";
	style.description = "q_{2}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62092:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the vertex [(q+1000)#times10]";
	style.description = "q_{3}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62093:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the vertex [(q+1000)#times10]";
	style.description = "q_{4}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62094:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 1 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck0}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 30, 0, 150000, "Np LFC"}, {31, 61, 0, 150000, "Np LPC"}, {62, 92, 0, 150000, "PorLNp LFC"}, {93, 123, 0, 150000, "PorLNp LPC"}, {124, 154, 0, 150000, "WCTrk0p LFC"}, {155, 185, 0, 150000, "WCTrk0p LPC"}, {186, 216, 0, 150000, "WCG0p LFC"}, {217, 247, 0, 150000, "WCG0p LPC"} };
	return true;
  case 62095:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 2 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck1}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 30, 0, 120000, "Np LFC"}, {31, 61, 0, 120000, "Np LPC"}, {62, 92, 0, 120000, "PorLNp LFC"}, {93, 123, 0, 120000, "PorLNp LPC"}, {124, 154, 0, 120000, "WCTrk0p LFC"}, {155, 185, 0, 120000, "WCTrk0p LPC"}, {186, 216, 0, 120000, "WCG0p LFC"}, {217, 247, 0, 120000, "WCG0p LPC"} };
	return true;
  case 62096:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 3 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck2}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62097:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 4 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck3}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62098:   // sp11_wcx0pNp_val
	style.title_axis_user = "Charge of muon spacepoint 5 from the muon end [(q+1000)#times10]";
	style.description = "q_{bck4}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 25, 0, 100000, "Np LFC"}, {26, 51, 0, 100000, "Np LPC"}, {52, 77, 0, 100000, "PorLNp LFC"}, {78, 103, 0, 100000, "PorLNp LPC"}, {104, 129, 0, 100000, "WCTrk0p LFC"}, {130, 155, 0, 100000, "WCTrk0p LPC"}, {156, 181, 0, 100000, "WCG0p LFC"}, {182, 207, 0, 100000, "WCG0p LPC"} };
	return true;
  case 62099:   // sp11_wcx0pNp_val
	style.title_axis_user = "Median muon spacepoint charge [(q+1000)#times10]";
	style.description = "q_{med}, all channels (sp11_wcx0pNp_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 28, 18000, 46000, "Np LFC"}, {29, 57, 18000, 46000, "Np LPC"}, {58, 86, 18000, 46000, "PorLNp LFC"}, {87, 115, 18000, 46000, "PorLNp LPC"}, {116, 144, 18000, 46000, "WCTrk0p LFC"}, {145, 173, 18000, 46000, "WCTrk0p LPC"}, {174, 202, 18000, 46000, "WCG0p LFC"}, {203, 231, 18000, 46000, "WCG0p LPC"} };
	return true;
  case 63001:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3} 0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 24, 0, 360000, ""} };
	return true;
  case 63002:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3} 0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 24, 0, 360000, ""} };
	return true;
  case 63003:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3} WCG0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 24, 0, 360000, ""} };
	return true;
  case 63004:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3} WCG0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 24, 0, 360000, ""} };
	return true;
  case 63005:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3} WCTrk0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 24, 0, 360000, ""} };
	return true;
  case 63006:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3} WCTrk0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 24, 0, 360000, ""} };
	return true;
  case 63007:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3} 0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 63008:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3} 0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 63009:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3} WCG0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 63010:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3} WCG0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 63011:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3} WCTrk0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 63012:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3} WCTrk0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 63013:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3} 0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 63014:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3} 0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 63015:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3} WCG0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 63016:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3} WCG0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 63017:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3} WCTrk0p LFC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 63018:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3} WCTrk0p LPC (sp0p_N3_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 63019:   // sp0p_N3_val
	style.title_axis_user = "Sum of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum3}, all channels (sp0p_N3_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 24, 0, 360000, "0p LFC"}, {25, 49, 0, 360000, "0p LPC"}, {50, 74, 0, 360000, "WCG0p LFC"}, {75, 99, 0, 360000, "WCG0p LPC"}, {100, 124, 0, 360000, "WCTrk0p LFC"}, {125, 149, 0, 360000, "WCTrk0p LPC"} };
	return true;
  case 63020:   // sp0p_N3_val
	style.title_axis_user = "Largest of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max3}, all channels (sp0p_N3_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 28, 0, 210000, "0p LFC"}, {29, 57, 0, 210000, "0p LPC"}, {58, 86, 0, 210000, "WCG0p LFC"}, {87, 115, 0, 210000, "WCG0p LPC"}, {116, 144, 0, 210000, "WCTrk0p LFC"}, {145, 173, 0, 210000, "WCTrk0p LPC"} };
	return true;
  case 63021:   // sp0p_N3_val
	style.title_axis_user = "Largest neighbour difference of the first 3 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max3}, all channels (sp0p_N3_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 30, 0, 150000, "0p LFC"}, {31, 61, 0, 150000, "0p LPC"}, {62, 92, 0, 150000, "WCG0p LFC"}, {93, 123, 0, 150000, "WCG0p LPC"}, {124, 154, 0, 150000, "WCTrk0p LFC"}, {155, 185, 0, 150000, "WCTrk0p LPC"} };
	return true;
  case 64001:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5} 0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 20, 30000, 430000, ""} };
	return true;
  case 64002:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5} 0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 20, 30000, 430000, ""} };
	return true;
  case 64003:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5} WCG0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 20, 30000, 430000, ""} };
	return true;
  case 64004:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5} WCG0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 20, 30000, 430000, ""} };
	return true;
  case 64005:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5} WCTrk0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 20, 30000, 430000, ""} };
	return true;
  case 64006:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5} WCTrk0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 20, 30000, 430000, ""} };
	return true;
  case 64007:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5} 0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 64008:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5} 0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 64009:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5} WCG0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 64010:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5} WCG0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 64011:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5} WCTrk0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 64012:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5} WCTrk0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 28, 0, 210000, ""} };
	return true;
  case 64013:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5} 0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 64014:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5} 0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 64015:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5} WCG0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 64016:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5} WCG0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 64017:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5} WCTrk0p LFC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 64018:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5} WCTrk0p LPC (sp0p_N5_val)";
	style.segment_divisions = 505;
	style.segments = { {0, 30, 0, 150000, ""} };
	return true;
  case 64019:   // sp0p_N5_val
	style.title_axis_user = "Sum of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{sum5}, all channels (sp0p_N5_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 20, 30000, 430000, "0p LFC"}, {21, 41, 30000, 430000, "0p LPC"}, {42, 62, 30000, 430000, "WCG0p LFC"}, {63, 83, 30000, 430000, "WCG0p LPC"}, {84, 104, 30000, 430000, "WCTrk0p LFC"}, {105, 125, 30000, 430000, "WCTrk0p LPC"} };
	return true;
  case 64020:   // sp0p_N5_val
	style.title_axis_user = "Largest of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "q_{max5}, all channels (sp0p_N5_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 28, 0, 210000, "0p LFC"}, {29, 57, 0, 210000, "0p LPC"}, {58, 86, 0, 210000, "WCG0p LFC"}, {87, 115, 0, 210000, "WCG0p LPC"}, {116, 144, 0, 210000, "WCTrk0p LFC"}, {145, 173, 0, 210000, "WCTrk0p LPC"} };
	return true;
  case 64021:   // sp0p_N5_val
	style.title_axis_user = "Largest neighbour difference of the first 5 muon spacepoint charges [(q+1000)#times10]";
	style.description = "dq_{max5}, all channels (sp0p_N5_val)";
	style.segment_divisions = 503;
	style.segments = { {0, 30, 0, 150000, "0p LFC"}, {31, 61, 0, 150000, "0p LPC"}, {62, 92, 0, 150000, "WCG0p LFC"}, {93, 123, 0, 150000, "WCG0p LPC"}, {124, 154, 0, 150000, "WCTrk0p LFC"}, {155, 185, 0, 150000, "WCTrk0p LPC"} };
	return true;
  default:
	return false;
  }
}
