/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2016 team free-astro (see more in AUTHORS file)
 * Reference site is http://free-astro.vinvin.tf/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/param.h>
#include <fcntl.h>
#if (defined(__APPLE__) && defined(__MACH__))
#include <mach/task.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <mach/vm_statistics.h>
#endif
#include <dirent.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#define _WITH_GETLINE
#include <stdio.h>
#include <string.h>
#include <locale.h>
#include <libconfig.h>
#include <assert.h>
#include <math.h>
#include "core/siril.h"
#include "core/proto.h"
#include "io/conversion.h"
#include "gui/callbacks.h"
#include "io/single_image.h"

static const char *keywords[] = { "working-directory", "libraw-settings",
		"debayer-settings", "registration-settings", "stacking-settings",
		"misc-settings" };
enum token_index {
	WD = 0, RAW = 1, BAY = 2, REG = 3, STK = 4, MISC = 5, NOTOK
};

int round_to_int(double x) {
	assert(x >= INT_MIN-0.5);
	assert(x <= INT_MAX+0.5);
	if (x >= 0.0)
		return (int) (x + 0.5);
	return (int) (x - 0.5);
}

WORD round_to_WORD(double x) {
	if (x <= 0.0)
		return (WORD) 0;
	if (x > USHRT_MAX_DOUBLE)
		return USHRT_MAX;
	return (WORD) (x + 0.5);
}

BYTE round_to_BYTE(double x) {
	if (x <= 0.0)
		return (BYTE) 0;
	if (x > UCHAR_MAX_DOUBLE)
		return UCHAR_MAX;
	return (BYTE) (x + 0.5);
}

BYTE conv_to_BYTE(double x) {
	if (x == 0.0)
		return (BYTE) 0;
	if (x == USHRT_MAX_DOUBLE)
		return UCHAR_MAX;
	x = ((x / USHRT_MAX_DOUBLE) * UCHAR_MAX_DOUBLE);
	return((BYTE) (x));
}

/* returns TRUE if fit has 3 layers */
gboolean isrgb(fits *fit) {
	return (fit->naxis == 3);
}

/* returns TRUE if str ends with ending, case insensitive */
gboolean ends_with(const char *str, const char *ending) {
	if (!str || str[0] == '\0')
		return FALSE;
	if (!ending || ending[0] == '\0')
		return TRUE;
	int ending_len = strlen(ending);
	int str_len = strlen(str);
	if (ending_len > str_len)
		return FALSE;
	return !strncasecmp(str + str_len - ending_len, ending, ending_len);
}

/* searches for an extension '.something' in filename from the end, and returns
 * the index of the first '.' found */
int get_extension_index(const char *filename) {
	int i;
	if (filename == NULL || filename[0] == '\0')
		return -1;
	i = strlen(filename) - 1;
	do {
		if (filename[i] == '.')
			return i;
		i--;
	} while (i > 0);
	return -1;
}

/* Get the extension of a file, without the dot.
 * The returned pointed is from the filename itself or NULL. */
const char *get_filename_ext(const char *filename) {
	const char *dot = strrchr(filename, '.');
	if (!dot || dot == filename)
		return NULL;
	return dot + 1;
}

// tests whether the given file is either regular or a symlink
// Return value is 1 if file is readable (not actually opened to verify)
int is_readable_file(const char *filename) {
	struct stat sts;
	if (stat(filename, &sts))
		return 0;
	if (S_ISREG (sts.st_mode) || S_ISLNK(sts.st_mode))
		return 1;
	return 0;
}

/* Tests if filename is the canonical name of a known file type
 * `type' is set according to the result of the test,
 * `realname' (optionnal) is set according to the found file name.
 * If filename contains an extension, only this file name is tested, else all
 * extensions are tested for the file name until one is found. */
int stat_file(const char *filename, image_type *type, char *realname) {
	// FIXME: realname should be a char ** to be allocated in the function
	char *test_name;
	int k;
	const char *ext = get_filename_ext(filename);
	*type = TYPEUNDEF;	// default value

	/* check for an extension in filename and isolate it, including the . */
	if (!filename || filename[0] == '\0')
		return 1;
	/* if filename has an extension, we only test for it */
	if (ext) {
		if (is_readable_file(filename)) {
			if (realname)
				strcpy(realname, filename);
			*type = get_type_for_extension(ext);
			return 0;
		}
		return 1;
	}

	test_name = malloc(strlen(filename) + 10);
	/* else, we can test various file extensions */
	/* first we test lowercase, then uppercase */
	for (k = 0; k < 2; k++) {
		int i = 0;
		while (supported_extensions[i]) {
			char *str;
			if (k == 0)
				str = strdup(supported_extensions[i]);
			else str = convtoupper(supported_extensions[i]);
			snprintf(test_name, 255, "%s%s", filename, str);
			free(str);
			if (is_readable_file(test_name)) {
				*type = get_type_for_extension(supported_extensions[i]+1);
				assert(*type != TYPEUNDEF);
				if (realname)
					strcpy(realname, test_name);
				return 0;
			}
			i++;
		}
	}
	return 1;
}

/******************************** CONFIGURATION FILE **********************************
 *                       We use libconf to generate this file                         *
 * ************************************************************************************/

int readinitfile() {
	config_t config;
	const char *dir = NULL, *swap_dir = NULL, *extension = NULL;

	if (!com.initfile)
		return 1;

	config_init(&config);

	if (config_read_file(&config, com.initfile) == CONFIG_FALSE)
		return 1;
	siril_log_message("Loading init file: %s\n", com.initfile);

	/* Working directory */
	if (config_lookup_string(&config, keywords[WD], &dir)) {
		if (changedir(dir)) {
			siril_log_message(
					"Reverting current working directory to startup directory, the saved directory is not available anymore\n");
			set_GUI_CWD();
			writeinitfile();
		}
	}

	/* Libraw setting */
	config_setting_t *raw_setting = config_lookup(&config, keywords[RAW]);
	if (raw_setting) {
		config_setting_lookup_float(raw_setting, "mul_0", &com.raw_set.mul[0]);
		config_setting_lookup_float(raw_setting, "mul_2", &com.raw_set.mul[2]);
		config_setting_lookup_float(raw_setting, "bright", &com.raw_set.bright);
		config_setting_lookup_int(raw_setting, "auto", &com.raw_set.auto_mul);
		config_setting_lookup_int(raw_setting, "cam_wb",
				&com.raw_set.use_camera_wb);
		config_setting_lookup_int(raw_setting, "auto_wb",
				&com.raw_set.use_auto_wb);
		config_setting_lookup_int(raw_setting, "user_qual",
				&com.raw_set.user_qual);
		config_setting_lookup_float(raw_setting, "gamm_0",
				&com.raw_set.gamm[0]);
		config_setting_lookup_float(raw_setting, "gamm_1",
				&com.raw_set.gamm[1]);
		config_setting_lookup_int(raw_setting, "user_black",
				&com.raw_set.user_black);
	}

	/* Debayer setting */
	config_setting_t *debayer_setting = config_lookup(&config, keywords[BAY]);
	if (debayer_setting) {
		config_setting_lookup_bool(debayer_setting, "ser_use_bayer_header",
				&com.debayer.ser_use_bayer_header);
		config_setting_lookup_int(debayer_setting, "pattern",
				&com.debayer.bayer_pattern);
		int inter;
		config_setting_lookup_int(debayer_setting, "inter", &inter);
		com.debayer.bayer_inter = inter;
	}

	/* Registration setting */
	config_setting_t *reg_setting = config_lookup(&config, keywords[REG]);
	if (reg_setting) {
		config_setting_lookup_int(reg_setting, "method", &com.reg_settings);
	}

	/* Stacking setting */
	config_setting_t *stack_setting = config_lookup(&config, keywords[STK]);
	if (stack_setting) {
		config_setting_lookup_int(stack_setting, "method", &com.stack.method);
		config_setting_lookup_int(stack_setting, "rejection",
				&com.stack.rej_method);
		config_setting_lookup_int(stack_setting, "normalisation",
				&com.stack.normalisation_method);
		config_setting_lookup_float(stack_setting, "maxmem",
				&com.stack.memory_percent);
	}
	if (com.stack.memory_percent <= 0.0001)
		com.stack.memory_percent = 0.9;

	/* Misc setting */
	config_setting_t *misc_setting = config_lookup(&config, keywords[MISC]);
	if (misc_setting) {
		config_setting_lookup_bool(misc_setting, "confirm", &com.dontShowConfirm);
		config_setting_lookup_string(misc_setting, "swap_directory", &swap_dir);
		config_setting_lookup_string(misc_setting, "extension", &extension);
		set_GUI_misc();
	}
	if (swap_dir && swap_dir[0] != '\0') {
		if (com.swap_dir)
			free(com.swap_dir);
		com.swap_dir = strdup(swap_dir);
	}
	else {
		const char* sw_dir = g_get_tmp_dir();
		com.swap_dir = strdup(sw_dir);
	}
	if (extension && extension[0] != '\0') {
		if (com.ext)
			free(com.ext);
		com.ext = strdup(extension);
	}
	else {
		com.ext = strdup(".fit");
	}
	return 0;
}

static void _save_wd(config_t *config, config_setting_t *root) {
	config_setting_t *directory;

	directory = config_setting_add(root, keywords[WD], CONFIG_TYPE_STRING);
	config_setting_set_string(directory, com.wd);
}

static void _save_libraw(config_t *config, config_setting_t *root) {
	config_setting_t *libraw_group, *raw_setting;

	libraw_group = config_setting_add(root, keywords[RAW], CONFIG_TYPE_GROUP);

	raw_setting = config_setting_add(libraw_group, "mul_0", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.mul[0]);

	raw_setting = config_setting_add(libraw_group, "mul_2", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.mul[2]);

	raw_setting = config_setting_add(libraw_group, "bright", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.bright);

	raw_setting = config_setting_add(libraw_group, "auto", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.auto_mul);

	raw_setting = config_setting_add(libraw_group, "cam_wb", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.use_camera_wb);

	raw_setting = config_setting_add(libraw_group, "auto_wb", CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.use_auto_wb);

	raw_setting = config_setting_add(libraw_group, "user_qual",
			CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.user_qual);

	raw_setting = config_setting_add(libraw_group, "gamm_0", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.gamm[0]);

	raw_setting = config_setting_add(libraw_group, "gamm_1", CONFIG_TYPE_FLOAT);
	config_setting_set_float(raw_setting, com.raw_set.gamm[1]);

	raw_setting = config_setting_add(libraw_group, "user_black",
			CONFIG_TYPE_INT);
	config_setting_set_int(raw_setting, com.raw_set.user_black);
}

static void _save_debayer(config_t *config, config_setting_t *root) {
	config_setting_t *debayer_group, *debayer_setting;

	debayer_group = config_setting_add(root, keywords[BAY], CONFIG_TYPE_GROUP);

	debayer_setting = config_setting_add(debayer_group, "ser_use_bayer_header",
			CONFIG_TYPE_BOOL);
	config_setting_set_bool(debayer_setting, com.debayer.ser_use_bayer_header);

	debayer_setting = config_setting_add(debayer_group, "pattern", CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.debayer.bayer_pattern);

	debayer_setting = config_setting_add(debayer_group, "inter", CONFIG_TYPE_INT);
	config_setting_set_int(debayer_setting, com.debayer.bayer_inter);
}

static void _save_registration(config_t *config, config_setting_t *root) {
	config_setting_t *reg_group, *reg_setting;

	reg_group = config_setting_add(root, keywords[REG], CONFIG_TYPE_GROUP);

	reg_setting = config_setting_add(reg_group, "method", CONFIG_TYPE_INT);
	config_setting_set_int(reg_setting, com.reg_settings);
}

static void _save_stacking(config_t *config, config_setting_t *root) {
	config_setting_t *stk_group, *stk_setting;

	stk_group = config_setting_add(root, keywords[STK], CONFIG_TYPE_GROUP);

	stk_setting = config_setting_add(stk_group, "method", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.stack.method);

	stk_setting = config_setting_add(stk_group, "rejection", CONFIG_TYPE_INT);
	config_setting_set_int(stk_setting, com.stack.rej_method);

	stk_setting = config_setting_add(stk_group, "maxmem", CONFIG_TYPE_FLOAT);
	config_setting_set_float(stk_setting, com.stack.memory_percent);
}

static void _save_misc(config_t *config, config_setting_t *root) {
	config_setting_t *misc_group, *misc_setting;

	misc_group = config_setting_add(root, keywords[MISC], CONFIG_TYPE_GROUP);

	misc_setting = config_setting_add(misc_group, "swap_directory", CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.swap_dir);

	misc_setting = config_setting_add(misc_group, "extension", CONFIG_TYPE_STRING);
	config_setting_set_string(misc_setting, com.ext);

	misc_setting = config_setting_add(misc_group, "confirm", CONFIG_TYPE_BOOL);
	config_setting_set_bool(misc_setting, com.dontShowConfirm);
}

int writeinitfile() {
	config_t config;
	config_setting_t *root;

	config_init(&config);
	root = config_root_setting(&config);

	_save_wd(&config, root);
	_save_libraw(&config, root);
	_save_debayer(&config, root);
	_save_registration(&config, root);
	_save_stacking(&config, root);
	_save_misc(&config, root);

	if (!config_write_file(&config, com.initfile)) {
		fprintf(stderr, "Error while writing file.\n");
		config_destroy(&config);
		return 1;
	}
	config_destroy(&config);

	return 0;
}

int checkinitfile() {
	char *home;
	struct stat sts;

	// try to read the file given on command line
	if (com.initfile && !readinitfile()) {
		return 0;
	}

	// no file given on command line, set initfile to default location
	if ((home = getenv("HOME")) == NULL) {
		fprintf(stderr,
				"Could not get the environment variable $HOME, no config file.\n");
		return 1;
	}
	com.initfile = malloc(strlen(home) + 20);
	sprintf(com.initfile, "%s/.siril/siril.cfg", home);
	if (readinitfile()) {	// couldn't read it
		char filename[255];

		set_GUI_CWD();
		// if that fails, check and create the default ini file
		snprintf(filename, 255, "%s/.siril", home);
		if (stat(filename, &sts) != 0) {
			if (errno == ENOENT) {
				if (mkdir(filename, 0755)) {
					fprintf(stderr, "Could not create dir %s, please check\n",
							filename);
					return 1;
				}
				com.swap_dir = strdup (g_get_tmp_dir());
				com.ext = strdup(".fit");
				return (writeinitfile());
			}
		}

		if (!(S_ISDIR(sts.st_mode))) {
			fprintf(stderr,
					"There is a file named %s, that is not a directory.\n"
					"Remove or rename it first\n",
					filename);
			return 1;
		}

		com.swap_dir = strdup (g_get_tmp_dir());
		com.ext = strdup(".fit");
		return (writeinitfile());
	}
	return 0;
}
/***************************************************************************************/

/* Try to change the CWD to the argument, absolute or relative.
 * If success, the new CWD is written to com.wd */
int changedir(const char *dir) {
	if (dir == NULL || dir[0] == '\0')
		return 1;
	if (!chdir(dir)) {
		char str[256];

		/* do we need to search for sequences in the directory now? We still need to
		 * press the check seq button to display the list, and this is also done there. */
		/* check_seq();
		 update_sequence_list();*/
		if (dir[0] == '/') {
			if (com.wd)
				free(com.wd);
			com.wd = strdup(dir);
			if (!com.wd)
				return 1;
		} else {
			// dir can be a relative path
			com.wd = realloc(com.wd, PATH_MAX);
			if (!com.wd)
				return 1;
			com.wd = getcwd(com.wd, PATH_MAX);
		}
		siril_log_message("Setting CWD (Current Working Directory) to `%s'\n",
				com.wd);
		set_GUI_CWD();

		snprintf(str, 255, "%s v%s - %s", PACKAGE, VERSION, dir);
		gtk_window_set_title(
				GTK_WINDOW(gtk_builder_get_object(builder, "main_window")),
				str);
		return 0;
	}
	siril_log_message("Could not change directory.\n");
	return 1;
}

/* This method populates the sequence combo box with the sequences found in the CWD.
 * If only one sequence is found, or if a sequence whose name matches the
 * possibly NULL argument is found, it is automatically selected, which triggers
 * its loading */
int update_sequences_list(const char *sequence_name_to_select) {
	GtkComboBoxText *seqcombo;
	DIR *dir;
	struct dirent *file;
	char filename[256];
	int number_of_loaded_sequences = 0;
	int index_of_seq_to_load = -1;

	if ((dir = opendir(com.wd)) == NULL) {
		fprintf(stderr, "working directory cannot be opened.\n");
		com.wd[0] = '\0';
		return 1;
	}

	// clear the previous list
	seqcombo = GTK_COMBO_BOX_TEXT(
			gtk_builder_get_object(builder, "sequence_list_combobox"));
	gtk_combo_box_text_remove_all(seqcombo);

	while ((file = readdir(dir)) != NULL) {
		char *suf;

		if ((suf = strstr(file->d_name, ".seq")) && strlen(suf) == 4) {
			sequence *seq = readseqfile(file->d_name);
			if (seq != NULL) {
				strncpy(filename, file->d_name, 255);
				free_sequence(seq, TRUE);
				gtk_combo_box_text_append_text(seqcombo, filename);
				if (sequence_name_to_select
						&& !strncmp(filename, sequence_name_to_select,
								strlen(filename))) {
					index_of_seq_to_load = number_of_loaded_sequences;
				}
				++number_of_loaded_sequences;
			}
		}
	}
	closedir(dir);
	if (!number_of_loaded_sequences) {
		fprintf(stderr, "No valid sequence found in CWD.\n");
		return -1;
	} else {
		fprintf(stdout, "Loaded %d sequence(s)\n", number_of_loaded_sequences);
	}

	if (number_of_loaded_sequences > 1 && index_of_seq_to_load < 0) {
		//g_signal_handlers_block_by_func(GTK_WIDGET(seqcombo), on_seqproc_entry_changed, NULL);
		gtk_combo_box_popup(GTK_COMBO_BOX(seqcombo));
		//g_signal_handlers_unblock_by_func(GTK_WIDGET(seqcombo), on_seqproc_entry_changed, NULL);
	} else if (index_of_seq_to_load >= 0)
		gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), index_of_seq_to_load);
	else
		gtk_combo_box_set_active(GTK_COMBO_BOX(seqcombo), 0);
	return 0;
}

void update_used_memory() {
#if defined(__linux__)
	unsigned long size, resident, share, text, lib, data, dt;
	static int page_size_in_k = 0;
	const char* statm_path = "/proc/self/statm";
	FILE *f = fopen(statm_path, "r");

	if (page_size_in_k == 0) {
		page_size_in_k = getpagesize() / 1024;
	}
	if (!f) {
		perror(statm_path);
		return;
	}
	if (7 != fscanf(f, "%lu %lu %lu %lu %lu %lu %lu",
			&size, &resident, &share, &text, &lib, &data, &dt)) {
		perror(statm_path);
		fclose(f);
		return;
	}
	fclose(f);
	set_GUI_MEM(resident * page_size_in_k);
#elif (defined(__APPLE__) && defined(__MACH__))
	struct task_basic_info t_info;

	mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
	set_GUI_MEM((unsigned long) t_info.resident_size / 1024UL);
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). In fact, it could work with linux */
	struct rusage usage;

	getrusage(RUSAGE_SELF, &usage);
	set_GUI_MEM((unsigned long) usage.ru_maxrss);
#else
	set_GUI_MEM((unsigned long) 0);
#endif
}

int get_available_memory_in_MB() {
	int mem = 2048; /* this is the default value if we can't retrieve any values */

#if defined(__linux__)
	FILE* fp = fopen("/proc/meminfo", "r");
	if (fp != NULL) {
		size_t bufsize = 1024 * sizeof(char);
		char* buf = (char*) malloc(bufsize);
		long value = -1L;
		while (getline(&buf, &bufsize, fp) >= 0) {
			if (strncmp(buf, "MemAvailable", 12) != 0)
				continue;
			sscanf(buf, "%*s%ld", &value);
			break;
		}
		fclose(fp);
		free((void*) buf);
		if (value != -1L)
			mem = (int) (value / 1024L);
	}
#elif (defined(__APPLE__) && defined(__MACH__))
	vm_size_t page_size;
	mach_port_t mach_port;
	mach_msg_type_number_t count;
	vm_statistics64_data_t vm_stats;

	mach_port = mach_host_self();
	count = sizeof(vm_stats) / sizeof(natural_t);
	if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
			KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
					(host_info64_t)&vm_stats, &count))	{

		int64_t unused_memory = ((int64_t)vm_stats.free_count +
				(int64_t)vm_stats.inactive_count +
				(int64_t)vm_stats.wire_count) * (int64_t)page_size;

		mem = (int) ((unused_memory) / (1024 * 1024));
	}
#elif defined(BSD) /* BSD (DragonFly BSD, FreeBSD, OpenBSD, NetBSD). ----------- */
	FILE* fp = fopen("/var/run/dmesg.boot", "r");
	if (fp != NULL) {
		size_t bufsize = 1024 * sizeof(char);
		char* buf = (char*) malloc(bufsize);
		long value = -1L;
		while (getline(&buf, &bufsize, fp) >= 0) {
			if (strncmp(buf, "avail memory", 12) != 0)
				continue;
			sscanf(buf, "%*s%*s%*s%ld", &value);
			break;
		}
		fclose(fp);
		free((void*) buf);
		if (value != -1L)
			mem = (int) (value / 1024L);
	}
#else

#endif
	return mem;
}

#if 0
/* returns true if the command theli is available */
gboolean theli_is_available() {
	int retval = system("theli > /dev/null");
	if (WIFEXITED(retval))
	return 0 == WEXITSTATUS(retval); // 0 if it's available, 127 if not
	return FALSE;
}
#endif

/* expands the ~ in filenames */
void expand_home_in_filename(char *filename, int size) {
	if (filename[0] == '~' && filename[1] == '\0')
		strcat(filename, "/");
	int len = strlen(filename);
	if (len < 2)
		return;		// not very necessary now with the first line
	if (filename[0] == '~' && filename[1] == '/') {
		char *homepath = getenv("HOME");
		int j, homelen = strlen(homepath);
		if (len + homelen > size - 1) {
			siril_log_message("file name is too long, not expanding it\n");
			return;
		}
		for (j = len; j > 0; j--)		// edit in place
			filename[j + homelen - 1] = filename[j];
		// the -1 above is tricky: it's the removal of the ~ character from
		// the original string
		strncpy(filename, homepath, homelen);
	}
}

WORD get_normalized_value(fits *fit) {
	image_find_minmax(fit, 0);
	if (fit->maxi <= UCHAR_MAX)
		return UCHAR_MAX;
	return USHRT_MAX;
}

/* This function reads a text file and displays it in the
 * show_data_dialog */
void read_and_show_textfile(char *path, char *title) {
	char line[64] = "";
	char txt[1024] = "";

	FILE *f = fopen(path, "r");
	if (!f) {
		show_dialog("File not found", "Error", "gtk-dialog-error");
		return;
	}
	while (fgets(line, sizeof(line), f) != NULL)
		strcat(txt, line);
	show_data_dialog(txt, title);
	fclose(f);
}

/* Exchange the two parameters of the function 
 * Usefull in Dynamic PSF (PSF.c) */
void swap_param(double *a, double *b) {
	double tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

// in-place quick sort, of array a of size n
void quicksort_d(double *a, int n) {
	if (n < 2)
		return;
	double p = a[n / 2];
	double *l = a;
	double *r = a + n - 1;
	while (l <= r) {
		if (*l < p) {
			l++;
			continue;
		}
		if (*r > p) {
			r--;
			continue; // we need to check the condition (l <= r) every time we change the value of l or r
		}
		double t = *l;
		*l++ = *r;
		*r-- = t;
	}
	quicksort_d(a, r - a + 1);
	quicksort_d(l, a + n - l);
}

// in-place quick sort, of array a of size n
void quicksort_s(WORD *a, int n) {
	if (n < 2)
		return;
	WORD p = a[n / 2];
	WORD *l = a;
	WORD *r = a + n - 1;
	while (l <= r) {
		if (*l < p) {
			l++;
			continue;
		}
		if (*r > p) {
			r--;
			continue; // we need to check the condition (l <= r) every time we change the value of l or r
		}
		WORD t = *l;
		*l++ = *r;
		*r-- = t;
	}
	quicksort_s(a, r - a + 1);
	quicksort_s(l, a + n - l);
}

double get_median_value_from_sorted_word_data(WORD *data, int size) {
	double median = 0;

	switch (size % 2) {
	case 0:
		median = 0.5 * (data[(size / 2) - 1] + data[size / 2]);
		break;
	default:
		median = data[(size - 1) / 2];
		break;
	}
	return median;
}

/* return the sigma of a set of data. In the same time it computes
 * the mean value */
double get_standard_deviation(WORD *data, int size) {
	int i;
	double sigma, mean;

	mean = sigma = 0.0;
	if (size < 2)
		return 0.0;
	for (i = 0; i < size; ++i) {
		mean += (double) data[i];
	}
	mean /= size;
	for (i = 0; i < size; ++i)
		sigma += pow((double) data[i] - mean, 2);
	sigma /= (size - 1);
	sigma = sqrt(sigma);
	return sigma;
}

/* returns a new string in uppercase */
char *convtoupper(char *str) {
	char *newstr, *p;
	p = newstr = strdup(str);
	while (*p) {
		*p = toupper(*p);
		++p;
	}

	return newstr;
}

/* This function returns a new string containing the dir name of a path.
 * Conversely, dirname modifies the argument.
 * The return value here must be freed. */
char *extract_path(const char *filename) {
	size_t pathlen;
	const char *p;
	char *dirname = NULL;

	p = strrchr(filename, '/');

	if (p == NULL) {
		dirname = malloc(1);
		dirname[0] = '\0';
		return dirname;
	}

	pathlen = p - filename;	// contains the trailing '/'
	dirname = malloc(pathlen + 1);
	strncpy(dirname, filename, pathlen);
	dirname[pathlen] = '\0';

	return dirname;
}

char *remove_ext_from_filename(const char *filename) {
	size_t filelen;
	const char *p;
	char *file = NULL;

	p = strrchr(filename, '.');

	if (p == NULL) {
		file = malloc(1);
		file[0] = '\0';
		return file;
	}

	filelen = p - filename;
	file = malloc(filelen + 1);
	strncpy(file, filename, filelen);
	file[filelen] = '\0';

	return file;
}

char *replace_spaces_from_filename(const char *filename) {
	GString *string = g_string_new(filename);
	size_t shift = 0;
	int i;

	for (i = 0; filename[i]; i++) {
		if (filename[i] == ' ') {
			string = g_string_insert_c(string, i + shift, '\\');
			shift++;
		}
	}
	return g_string_free(string, FALSE);
}

// append a string to the end of an existing string
char* str_append(char** data, const char* newdata) {
	char* p;
	int len = (*data ? strlen(*data) : 0);
	if ((p = realloc(*data, len + strlen(newdata) + 1)) == NULL) {
		free(p);
		siril_log_message("str_append: error allocating data\n");
		return NULL;
	}
	*data = p;
	strcpy(*data + len, newdata);
	return *data;
}

/* cut a base name to 120 characters and add a trailing underscore if needed.
 * WARNING: may return a newly allocated string and free the argument */
char *format_basename(char *root) {
	int len = strlen(root);
	if (len > 120) {
		root[120] = '\0';
		len = 120;
	}
	if (root[len-1] == '-' || root[len-1] == '_') {
		return root;
	}

	char *appended = malloc(len+2);
	sprintf(appended, "%s_", root);
	free(root);
	return appended;
}
