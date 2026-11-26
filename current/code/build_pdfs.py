"""Export independent PDF panels and compose the manuscript plates from those PDFs."""

import argparse
import json
import os
from itertools import groupby
from pathlib import Path
from shutil import copyfile
from tempfile import NamedTemporaryFile

import pymupdf as pdf


ROOT = Path(__file__).resolve().parents[1]
MASTERS = Path(os.environ.get('GPUV_FIGURE_MASTERS', ROOT / 'data' / 'masters'))


def panel_layout(key, width, height):
    """Coordinates in PDF points from the upper left. Each tile is assigned once."""
    layouts = {
        '1': {
            'a': [(0, 0, width, 282)], 'b': [(0, 282, 248, 479)],
            'c': [(248, 282, width, 479)],
            'd': [(0, 479, 248, 549), (0, 549, 343, 592), (0, 592, 353, height)],
            'e': [(248, 479, width, 549), (343, 549, width, 592), (353, 592, width, height)],
        },
        '2': {
            'a': [(0, 0, 237, 225)], 'b': [(237, 0, 425, 225)],
            'c': [(425, 0, width, 225)], 'd': [(0, 225, width, 437)],
            'e': [(0, 437, 152, height)], 'f': [(152, 437, 423, height)],
            'g': [(423, 437, width, height)],
        },
        '3': {'': [(0, 0, width, height)]},
        '4': {'a': [(0, 0, width, 283)], 'b': [(0, 283, 329, height)], 'c': [(329, 283, width, height)]},
        '5': {
            'a': [(0, 0, 300, 207)], 'b': [(300, 0, width, 207)],
            'c': [(0, 207, 300, 398)], 'd': [(300, 207, width, 398)],
            'e': [(0, 398, width, 529)], 'f': [(0, 529, width, height)],
        },
        '6': {
            'a': [(0, 0, 216, 215)], 'b': [(216, 0, 426, 215)],
            'c': [(426, 0, width, 215)], 'd': [(0, 215, 345, 440)],
            'e': [(0, 440, width, height)], 'f': [(345, 215, 508, 440)],
            'g': [(508, 215, width, 440)],
        },
        'S1': {
            'a': [(0, 0, 369, height)], 'b': [(369, 0, width, 153)],
            'c': [(369, 153, width, 345)], 'd': [(369, 345, width, height)],
        },
        'S2': {'a': [(0, 0, 449, 245)], 'b': [(449, 0, width, 245)], 'c': [(0, 245, width, height)]},
        'S3': {
            'a': [(0, 0, 239, 224)], 'b': [(239, 0, 458, 224)],
            'c': [(458, 0, width, 224)], 'd': [(0, 224, 234, 414)],
            'e': [(234, 224, width, 414)],
            'f': [(0, 414, 229, 499), (0, 499, 282, height)],
            'g': [(229, 414, 495, 499), (282, 499, 495, height)],
            'h': [(495, 414, width, height)],
        },
        'S4': {'a': [(0, 0, 258.5, height)], 'b': [(258.5, 0, width, height)]},
        'S5': {'a': [(0, 0, 188, height)], 'b': [(188, 0, 437, height)], 'c': [(437, 0, width, height)]},
        'S6': {
            'a': [(0, 0, 344, 375)], 'b': [(344, 34, 510, 480)],
            'c': [(510, 34, width, 480)], 'd': [(0, 375, 344, 480)],
            'legend': [(344, 0, width, 34)],
            'e': [(0, 480, 228, 728)], 'f': [(228, 480, 454, 728)],
            'g': [(454, 480, width, 728)], 'note': [(0, 728, width, height)],
        },
        'S7': {
            'a': [(0, 0, width, 177)], 'b': [(0, 177, 222, 337)],
            'c': [(222, 177, 368, 337)], 'd': [(368, 177, width, 337)],
            'e': [(0, 337, 262, 525)], 'f': [(262, 337, width, 525)],
            'g': [(0, 525, 333, height)], 'h': [(333, 525, width, height)],
        },
    }
    return {panel: [pdf.Rect(*box) for box in boxes] for panel, boxes in layouts[key].items()}


def vectorize_images(doc):
    """Turn existing chart pixels into independent filled paths in the same PDF layer."""
    for image in doc[0].get_image_info(xrefs=True):
        xref = image['xref']
        if not xref:
            raise ValueError('Unexpected inline image')
        pix = pdf.Pixmap(doc, xref)
        if pix.colorspace != pdf.csRGB:
            pix = pdf.Pixmap(pdf.csRGB, pix)
        if pix.alpha:
            raise ValueError('Unexpected alpha in chart image')
        width, height = pix.width, pix.height
        samples = pix.samples
        rows = [samples[y * pix.stride:(y + 1) * pix.stride] for y in range(height)]
        if width > 1000 or width == 932:
            if all(row == rows[0] for row in rows):
                rows = rows[:1]
        commands = []
        for y, row in enumerate(rows):
            pixels = [row[x:x + 3] for x in range(0, 3 * width, 3)]
            # Heatmap cells stay independent. Merge repeated colors only in color scales.
            runs = [(p, 1) for p in pixels] if (width, height) == (79, 17) else [
                (p, sum(1 for _ in run)) for p, run in groupby(pixels)
            ]
            x = 0
            for rgb, count in runs:
                r, g, b = (value / 255 for value in rgb)
                commands.append(f'{r:.7f} {g:.7f} {b:.7f} rg {x} {len(rows)-y-1} {count} 1 re f')
                x += count
        stream = f'q {1/width:.12f} 0 0 {1/len(rows):.12f} 0 0 cm\n' + '\n'.join(commands) + '\nQ\n'
        doc.update_object(xref, '<< /Type /XObject /Subtype /Form /BBox [0 0 1 1] /Resources << >> >>')
        doc.update_stream(xref, stream.encode('ascii'))


def save(doc, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=target.parent, suffix='.pdf', delete=False) as stream:
        pending = Path(stream.name)
    try:
        doc.save(pending, garbage=2, deflate=True)
        pending.replace(target)
    finally:
        pending.unlink(missing_ok=True)


def crop_panel(source, boxes, target, isolate=True):
    bounds = pdf.Rect(boxes[0])
    for box in boxes[1:]:
        bounds |= box
    origin = pdf.Point(bounds.x0, bounds.y0)
    output = pdf.open()
    page = output.new_page(width=bounds.width, height=bounds.height)
    for box in boxes:
        destination = box - pdf.Rect(origin, origin)
        if isolate:
            with pdf.open() as isolated:
                isolated.insert_pdf(source)
                isolated[0].clip_to_rect(box)
                page.show_pdf_page(destination, isolated, 0, clip=box)
        else:
            page.show_pdf_page(destination, source, 0, clip=box)
    save(output, target)
    output.close()
    return bounds


def standalone_city_map(target):
    """Restore the complete Fig. S3g map where the source plate overlaps S3f."""
    bounds = pdf.Rect(229, 414, 508, 668.02069)
    offset = pdf.Rect(bounds.x0, bounds.y0, bounds.x0, bounds.y0)
    with pdf.open() as doc, pdf.open(MASTERS / 'Fig_S3g_map_source.pdf') as map_doc:
        page = doc.new_page(width=bounds.width, height=bounds.height)
        clip = pdf.Rect(51, 87, 700, 654)
        scale = .4238
        region = pdf.Rect(208.37 + clip.x0*scale, 377.65 + clip.y0*scale,
                          208.37 + clip.x1*scale, 377.65 + clip.y1*scale)
        page.show_pdf_page(region-offset, map_doc, 0, clip=clip)
        font = Path(os.environ.get('ARIAL_BOLD_FILE', 'C:/Windows/Fonts/arialbd.ttf'))
        if not font.is_file():
            raise FileNotFoundError('Set ARIAL_BOLD_FILE to the installed Arial Bold font.')
        page.insert_font(fontname='ArialBold', fontfile=str(font))
        page.insert_text((3.5, 10.9), 'g', fontname='ArialBold', fontsize=8)
        doc.subset_fonts()
        save(doc, target)
    return bounds


def place_map_tiles(plate, panel, boxes, bounds):
    offset = pdf.Rect(bounds.x0, bounds.y0, bounds.x0, bounds.y0)
    for box in boxes:
        plate.show_pdf_page(box, panel, 0, clip=box-offset)


def assemble_s3_reference(source, path):
    """Share one resource set for the dense S3 panels in the reference assembly."""
    plate = source[0]
    size = plate.rect
    layout = panel_layout('S3', size.width, size.height)
    # Use the same coordinate system and font/form resources as the master.
    # Each region invokes the source content behind its own clipping path.
    plate.wrap_contents()
    content = plate.read_contents()
    page_to_pdf = ~plate.transformation_matrix
    tiled = []
    for boxes in layout.values():
        for box in boxes:
            area = box * page_to_pdf
            clip = f'q {area.x0:.7f} {area.y0:.7f} {area.width:.7f} {area.height:.7f} re W n\n'
            tiled.extend((clip.encode('ascii'), content, b'\nQ\n'))
    stream = source.get_new_xref()
    source.update_object(stream, '<<>>')
    source.update_stream(stream, b''.join(tiled))
    plate.set_contents(stream)
    save(source, path)


def describe(path):
    with pdf.open(path) as doc:
        page = doc[0]
        fonts = {name: bool(doc.extract_font(xref)[3]) for xref, _, _, name, *_ in page.get_fonts()}
        spans = [span for block in page.get_text('dict', flags=pdf.TEXTFLAGS_TEXT)['blocks']
                 for line in block.get('lines', []) for span in line['spans']]
        text_fonts = sorted({span['font'] for span in spans if span['text'].strip()})
        if not text_fonts or any('Arial' not in name for name in text_fonts):
            raise ValueError(f'{path.name}: unexpected text fonts {text_fonts}')
        image_count = len(page.get_image_info())
        if image_count:
            raise ValueError(f'{path.name}: found {image_count} raster images')
        return {'fonts': fonts, 'text_fonts': text_fonts, 'text_spans': len(spans),
                'images': image_count, 'size_pt': [page.rect.width, page.rect.height]}


def subtract_rectangle(area, hole):
    overlap = area & hole
    if overlap.is_empty:
        return [area]
    x0, y0, x1, y1 = area
    a, b, c, d = overlap
    return [box for box in (pdf.Rect(x0,y0,x1,b), pdf.Rect(x0,d,x1,y1),
                            pdf.Rect(x0,b,a,d), pdf.Rect(c,b,x1,d)) if not box.is_empty]


def prepare_master(key, source):
    if source[0].get_images(full=True):
        vectorize_images(source)
    if key not in ('5', 'S7'):
        return source

    # The revised CNS PDFs converted many unchanged labels to paths. Reuse
    # their editable versions while retaining the regions revised in CNS v3.
    editable = MASTERS / f'Fig_{key}_live_text_source.pdf'
    with pdf.open(editable) as previous:
        if previous[0].get_images(full=True):
            vectorize_images(previous)
        width, height = source[0].rect.width, source[0].rect.height
        result = pdf.open()
        page = result.new_page(width=width, height=height)
        if key == '5':
            revised = [pdf.Rect(300,207,width,398), pdf.Rect(573,137,width,159),
                       pdf.Rect(109.5,565.85,565,573.2), pdf.Rect(109.5,708.4,168,720.5)]
        else:
            revised = [pdf.Rect(0,521,width,height)]
        unchanged = [pdf.Rect(0,0,width,height)]
        for region in revised:
            unchanged = [part for box in unchanged for part in subtract_rectangle(box, region)]
        for region in unchanged:
            page.show_pdf_page(region, previous, 0, clip=region)
        for region in revised:
            page.show_pdf_page(region, source, 0, clip=region)
        return result


def build(key, output_dir, reuse_panels=False):
    audit = {}
    panel_paths = []
    with pdf.open(MASTERS / f'Fig_{key}_CNS.pdf') as master:
        source = master if reuse_panels else prepare_master(key, master)
        page_size = source[0].rect
        layout = panel_layout(key, page_size.width, page_size.height)
        plate_path = output_dir / 'figures' / f'Fig_{key}.pdf'
        if key == '3':
            path = output_dir / 'panels' / 'Fig_3' / 'Fig_3.pdf'
            if not reuse_panels:
                save(source, path)
            plate_path.parent.mkdir(parents=True, exist_ok=True)
            copyfile(path, plate_path)
            return plate_path, [path], {'Fig_3': describe(path), 'assembled': describe(plate_path)}
        with pdf.open() as combined:
            plate = combined.new_page(width=page_size.width, height=page_size.height)
            for panel, boxes in layout.items():
                layout_asset = panel in ('note', 'legend')
                name = f'Fig_{key}{"_" + panel if layout_asset else panel}'
                path = output_dir / ('layout_assets' if layout_asset else 'panels') / f'Fig_{key}' / f'{name}.pdf'
                if key == 'S3' and panel == 'g':
                    bounds = (pdf.Rect(229, 414, 508, page_size.height) if reuse_panels
                              else standalone_city_map(path))
                elif reuse_panels:
                    bounds = pdf.Rect(boxes[0])
                    for box in boxes[1:]:
                        bounds |= box
                else:
                    isolate = key not in ('5', 'S7') and (key, panel) != ('S3', 'h')
                    bounds = crop_panel(source, boxes, path, isolate=isolate)
                if key != 'S3' or reuse_panels:
                    with pdf.open(path) as independent_panel:
                        if key == 'S3' and panel == 'g':
                            place_map_tiles(plate, independent_panel, boxes, bounds)
                        else:
                            plate.show_pdf_page(bounds, independent_panel, 0)
                panel_paths.append(path)
                audit[name] = describe(path)
            if key == 'S3' and not reuse_panels:
                assemble_s3_reference(source, plate_path)
            else:
                save(combined, plate_path)
        if source is not master:
            source.close()

    audit['assembled'] = describe(plate_path)
    if key == 'S6':
        with pdf.open(MASTERS / 'Fig_S6_tree_key.pdf') as key_source, pdf.open() as legend:
            legend.insert_pdf(key_source, from_page=1, to_page=1)
            key_path = output_dir / 'layout_assets' / 'Fig_S6' / 'Fig_S6a_genus_key.pdf'
            save(legend, key_path)
        audit['genus_key'] = describe(key_path)
    return plate_path, panel_paths, audit


def main():
    keys = [str(i) for i in range(1, 7)] + [f'S{i}' for i in range(1, 8)]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fig', nargs='+', choices=keys, default=keys)
    parser.add_argument('--out', type=Path, default=ROOT / 'pdf')
    parser.add_argument('--reuse-panels', action='store_true', help='Only assemble existing panel PDFs')
    parser.add_argument('--collection-only', action='store_true', help='Bind the existing full figure PDFs')
    args = parser.parse_args()
    audit = {}
    with pdf.open() as all_figures:
        for key in args.fig:
            if args.collection_only:
                path = args.out / 'figures' / f'Fig_{key}.pdf'
            else:
                path, panels, audit[key] = build(key, args.out, args.reuse_panels)
                print(f'Fig {key}: {len(panels)} PDFs for assembly', flush=True)
            with pdf.open(path) as figure:
                all_figures.insert_pdf(figure)
        all_figures.set_toc([[1, f'Figure {key}', i + 1] for i, key in enumerate(args.fig)])
        name = 'All_figures' if args.fig == keys else 'Selected_figures'
        save(all_figures, args.out / f'{name}.pdf')
    if audit:
        (args.out / 'audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
