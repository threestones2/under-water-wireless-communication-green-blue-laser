"""Fill template via direct XML manipulation."""
import zipfile
import os
import shutil
from lxml import etree

TEMPLATE = r"D:\509所-余重俊（公开）\8.蓝绿课题\小论文\会议论文\conference-template-a4.docx"
OUTPUT   = r"D:\509所-余重俊（公开）\8.蓝绿课题\小论文\会议论文\APCCAS_SA-MCS_Paper.docx"
TEMP_DIR = r"D:\509所-余重俊（公开）\8.蓝绿课题\小论文\会议论文\temp_xml"

EM_DASH = "—"

TITLE = ("A Computationally Efficient Semi-Analytical Monte Carlo Method "
         "for Real-Time Underwater Optical Channel Estimation")

ABSTRACT_BODY = (
    "Real-time channel estimation is critical for adaptive underwater wireless "
    "optical communication (UWOC) systems deployed on autonomous underwater "
    "vehicles (AUVs), yet conventional Monte Carlo (MC) methods impose prohibitive "
    "computational costs. This paper presents a Semi-Analytical Monte Carlo "
    "Simulation (SA-MCS) framework that replaces physical photon reception with "
    "closed-form local estimation at each scattering node, eliminating the need "
    "for ray-marching and geometric hit-testing. Two scattering-angle sampling "
    "strategies are investigated: a look-up table (LUT) approach achieving O(1) "
    "constant-time sampling, and an inverse-transform method offering greater "
    "flexibility. The proposed methods are benchmarked against wavefront-coupled "
    "importance-sampling MC (WCI-MC) and ground-truth physical MC across three "
    "water types including clear ocean, coastal ocean, and turbid harbor, at "
    "photon counts from 10^3 to 10^5. Experimental results demonstrate that "
    "SA-MCS (LUT) achieves path loss estimation accuracy within 0.5 dB of the "
    "reference while reducing computation time by over an order of magnitude "
    "compared to WCI-MC. The LUT-based architecture, with its fixed memory "
    "footprint and constant-time sampling, is particularly well-suited for "
    "hardware acceleration on resource-constrained embedded platforms, enabling "
    "real-time channel adaptation for next-generation intelligent underwater "
    "circuits and systems."
)

KEYWORDS_BODY = (
    "underwater wireless optical communication, Monte Carlo simulation, "
    "semi-analytical method, channel estimation, edge computing, embedded systems"
)


def get_ns(root):
    """Get the w namespace from root."""
    for prefix, uri in root.nsmap.items():
        if 'wordprocessingml' in uri:
            return uri
    return 'http://purl.oclc.org/ooxml/wordprocessingml/main'


def extract_text(para, w_ns):
    """Get full text from a paragraph element."""
    texts = para.findall('.//{%s}t' % w_ns)
    return ''.join(t.text or '' for t in texts)


def set_para_text(para, text, w_ns):
    """Replace all text in a paragraph with new text, keeping first run's formatting."""
    runs = para.findall('{%s}r' % w_ns)

    if not runs:
        new_run = etree.SubElement(para, '{%s}r' % w_ns)
        new_t = etree.SubElement(new_run, '{%s}t' % w_ns)
        new_t.text = text
        new_t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
        return

    # Remove all runs except the first
    for run in runs[1:]:
        para.remove(run)

    # Replace text in first run
    first_run = runs[0]
    t_elements = first_run.findall('{%s}t' % w_ns)

    if t_elements:
        t_elements[0].text = text
        t_elements[0].set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
        for t in t_elements[1:]:
            first_run.remove(t)
    else:
        new_t = etree.SubElement(first_run, '{%s}t' % w_ns)
        new_t.text = text
        new_t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')


def main():
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)

    with zipfile.ZipFile(TEMPLATE, 'r') as zf:
        zf.extractall(TEMP_DIR)

    doc_path = os.path.join(TEMP_DIR, 'word', 'document.xml')
    tree = etree.parse(doc_path)
    root = tree.getroot()
    w_ns = get_ns(root)
    print(f"Using namespace: {w_ns}")

    body = root.find('{%s}body' % w_ns)
    paragraphs = body.findall('.//{%s}p' % w_ns)

    title_replaced = False
    abstract_replaced = False
    keywords_replaced = False

    for para in paragraphs:
        text = extract_text(para, w_ns).strip()

        if not title_replaced and "Paper Title" in text:
            set_para_text(para, TITLE, w_ns)
            title_replaced = True
            print("Title replaced.")

        if not abstract_replaced and text.startswith(f"Abstract{EM_DASH}"):
            new_text = f"Abstract{EM_DASH}{ABSTRACT_BODY}"
            set_para_text(para, new_text, w_ns)
            abstract_replaced = True
            print(f"Abstract replaced.")

        if not keywords_replaced and text.startswith(f"Keywords{EM_DASH}"):
            new_text = f"Keywords{EM_DASH}{KEYWORDS_BODY}"
            set_para_text(para, new_text, w_ns)
            keywords_replaced = True
            print(f"Keywords replaced.")

    if not title_replaced:
        print("ERROR: Title not found!")
    if not abstract_replaced:
        print("ERROR: Abstract not found!")
    if not keywords_replaced:
        print("ERROR: Keywords not found!")

    tree.write(doc_path, xml_declaration=True, encoding='UTF-8', standalone=True)

    if os.path.exists(OUTPUT):
        os.remove(OUTPUT)

    with zipfile.ZipFile(OUTPUT, 'w', zipfile.ZIP_DEFLATED) as zf_out:
        for root_dir, dirs, files in os.walk(TEMP_DIR):
            for f in files:
                file_path = os.path.join(root_dir, f)
                arcname = os.path.relpath(file_path, TEMP_DIR)
                zf_out.write(file_path, arcname)

    print(f"Saved to: {OUTPUT}")
    shutil.rmtree(TEMP_DIR)


if __name__ == "__main__":
    main()
