#ifndef CALLBACKS_H
#define CALLBACKS_H

#include <gtk/gtk.h>
#include <sys/time.h>
#include "core/siril.h"	// for sequence

void load_css_style_sheet (char *original, char *path);
void initialize_shortcuts();
void fill_about_dialog();
void refresh_GUI();
void initialize_remap();
char* siril_log_internal(const char* format, const char* color, va_list arglist);
char* siril_log_message(const char* format, ...);
char* siril_log_color_message(const char* format, const char* color, ...);
void show_time(struct timeval, struct timeval);
void set_cutoff_sliders_max_values();		// was set_upper_minmax
void set_cutoff_sliders_values();		// was set_ranges
void set_sliders_value_to_gfit();
void initialize_display_mode();
void set_display_mode();
void adjust_exclude(int n, gboolean changed);
void adjust_refimage(int n);
int adjust_sellabel();
void set_GUI_CWD();
void set_GUI_MEM(unsigned long size);
void test_and_allocate_reference_image(int vport);
void enable_view_reference_checkbox(gboolean status);
gboolean redraw(int vport, int remap);
void sliders_mode_set_state(sliders_mode);
void on_radiobutton_minmax_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_radiobutton_hilo_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void on_radiobutton_user_toggled(GtkToggleButton *togglebutton, gpointer user_data);
int copy_rendering_settings_when_chained(gboolean from_GUI);
void on_seqproc_entry_changed (GtkComboBox *widget,	gpointer user_data);

void update_libraw_interface();
void update_main_interface();
void set_libraw_settings_menu_available(gboolean);
void set_debayer_settings_menu_available(gboolean);
void clear_sampling_setting_box();
void set_GUI_CAMERA();
void set_GUI_LIBRAW();
void on_checkbutton_cam_toggled(GtkButton *button, gpointer user_data);
void on_checkbutton_auto_toggled(GtkButton *button, gpointer user_data);

typedef void (*selection_update_callback)();
void register_selection_update_callback(selection_update_callback f);
void unregister_selection_update_callback(selection_update_callback f);
void delete_selected_area();

int match_drawing_area_widget(GtkWidget *drawing_area, gboolean allow_rgb);
char *vport_number_to_name(int);
void calculate_fwhm(GtkWidget *);
void on_max_entry_changed(GtkEditable *editable, gpointer user_data);
void on_min_entry_changed(GtkEditable *editable, gpointer user_data);
void on_excludebutton_toggled(GtkToggleButton *togglebutton, gpointer user_data);
void display_filename();
void set_layers_for_assign();
void set_layers_for_registration();
void display_image_number(int index);
void on_ref_frame_toggled(GtkToggleButton *togglebutton, gpointer user_data);
gboolean on_drawingarea_key_press_event(GtkWidget *widget, GdkEventKey *event, gpointer user_data);
gboolean end_sequence_prepro(gpointer p);
void show_dialog(const char *text, const char *title, const char *icon);
void show_data_dialog(char *text, char *title);
void show_main_gray_window();
void show_rgb_window();
void hide_gray_window();
void hide_rgb_window();
int is_histogram_visible();
void set_cursor_waiting(gboolean waiting);

void start_in_new_thread(gpointer(*f)(gpointer p), gpointer p);
void stop_processing_thread();
void set_thread_run(gboolean b);
gboolean get_thread_run();
void progress_bar_reset_ready();
void set_progress_bar_data(const char *text, double percent);
void on_combozoom_changed(GtkComboBox *widget, gpointer user_data);
double get_zoom_val();
void zoomcombo_update_display_for_zoom();
void adjust_vport_size_to_image();
void scrollbars_hadjustment_changed_handler(GtkAdjustment *adjustment, gpointer user_data);
void scrollbars_vadjustment_changed_handler(GtkAdjustment *adjustment, gpointer user_data);
void set_output_filename_to_sequence_name();
void do_popup_rgbmenu(GtkWidget *my_widget, GdkEventButton *event);
void do_popup_graymenu(GtkWidget *my_widget, GdkEventButton *event);
void close_tab();
void activate_tab(int vport);
void control_window_switch_to_tab(main_tabs tab);
GtkWidget* lookup_widget (const gchar *widget_name);
void on_combodisplay_changed (GtkComboBox *widget, gpointer user_data);
void on_checkchain_toggled(GtkToggleButton *togglebutton, gpointer user_data);

void gtk_filter_add(GtkFileChooser *file_chooser, const gchar *title, const gchar *pattern, gboolean set_default);
void set_filters_dialog(GtkFileChooser *);
void set_prepro_button_sensitiveness();

void on_removegreen_activate(GtkMenuItem *menuitem, gpointer user_data);
void on_SCNR_Apply_clicked(GtkButton *button, gpointer user_data);
void fill_convert_list(GSList *list);

void update_spinCPU(int max);

#endif
