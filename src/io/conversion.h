#ifndef CONVERSION_H
#define CONVERSION_H

typedef struct {
	char *extension;			// name of the extension of raw
	char *manufacturer;			// name of the manufacturer
	sensor_pattern suggested_pattern;// type of bayer pattern. Not used for now
} supported_raw_list;

extern supported_raw_list supported_raw[];	//supported raw extensions
extern int get_nb_raw_supported();

void initialize_combocamera();
void list_format_available();
image_type get_type_for_extension_name(const char *extension);
void initialize_default_extension();
void check_for_conversion_form_completeness();
void initialize_converters();
int tofits(char *source, char *dest);
void update_raw_cfa_tooltip();
void update_statusbar_convert();
int count_selected_files();

#endif
