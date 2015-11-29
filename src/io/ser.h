#ifndef _SER_H_
#define _SER_H_

#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef enum {
	SER_MONO = 0,
	SER_BAYER_RGGB = 8,
	SER_BAYER_GRBG = 9,
	SER_BAYER_GBRG = 10,
	SER_BAYER_BGGR = 11,
	SER_BAYER_CYYM = 16,
	SER_BAYER_YCMY = 17,
	SER_BAYER_YMCY = 18,
	SER_BAYER_MYYC = 19,
	RGB = 100,	// SER v3
	BGR = 101	// SER v3
} ser_color;

typedef enum {
	SER_BIG_ENDIAN, SER_LITTLE_ENDIAN
} ser_endian;

typedef enum {
	SER_PIXEL_DEPTH_8 = 1, SER_PIXEL_DEPTH_16 = 2
} ser_pixdepth;

struct ser_struct {
	char *file_id;			// 14 byte
	int lu_id;			// 4	(14)
	ser_color color_id;		// 4	(18)
	ser_endian endianness;		// 4	(22)
	int image_width;		// 4	(26)
	int image_height;		// 4	(30)
	ser_pixdepth pixel_depth;	// 4	(34)
	unsigned int frame_count;	// 4	(38)
	char *observer;			// 40	(42)
	char *instrument;		// 40	(82)
	char *telescope;		// 40	(122)
	uint64_t date, date_utc;	// 8 and 8

	// internal stuff
	int fd;
	char *filename;
	unsigned int number_of_planes;
#ifdef _OPENMP
	omp_lock_t fd_lock;
#endif
};

void ser_init_struct(struct ser_struct *ser_file);
void ser_display_info(struct ser_struct *ser_file);
int ser_open_file(char *filename, struct ser_struct *ser_file);
int ser_write_header(struct ser_struct *ser_file);
int ser_close_file(struct ser_struct *ser_file);
int ser_read_frame(struct ser_struct *ser_file, int frame_no, fits *fit);
int ser_read_opened_partial(struct ser_struct *ser_file, int layer,
		int frame_no, WORD *buffer, const rectangle *area);
int ser_write_frame_from_fit(struct ser_struct *ser_file, fits *fit);
void set_combo_box_bayer_pattern(ser_color pattern);
void ser_manage_endianess_and_depth(struct ser_struct *ser_file, WORD *data, int frame_size);
void ser_manage_rgb_v3(struct ser_struct *ser_file, WORD *data, int frame_size);

#endif

