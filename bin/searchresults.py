import xmltodict


def print_hsp(hsp, query_id, subject_id, coverage, qlen):
    print(
        "\t".join([
            query_id,
            subject_id,
            'N/A',
            str(round(100 * coverage / qlen)),
            str(round(100 * (int(hsp['Hsp_query-to']) -
                             int(hsp['Hsp_query-from'])) / qlen))
        ]))


def tabulardatafromxml(inf):
    inf = open(inf)
    r = xmltodict.parse(inf.read(), attr_prefix='')
    print("# Fields: query id, subject id, subject sci names,"
          " % query coverage per subject, % query coverage per hsp")
    qlen = int(r['BlastOutput']['BlastOutput_query-len'])
    itrs = r['BlastOutput']['BlastOutput_iterations']
    if isinstance(itrs['Iteration'], list):
        for itr in itrs['Iteration']:
            process_iteration(itr, qlen)
    else:
        process_iteration(itrs['Iteration'], qlen)


def process_iteration(itr, qlen):
    query_id = itr['Iteration_query-def'].split(' ')[0]
    if isinstance(itr['Iteration_hits']['Hit'], list):
        print("# %d similar sequences found" % len(itr['Iteration_hits']['Hit']))
        for sseqs in itr['Iteration_hits']['Hit']:
            process_similar_sequence(sseqs, query_id, qlen)
    else:
        print("# 1 similar sequence found")
        process_similar_sequence(itr['Iteration_hits']['Hit'], query_id, qlen)

def process_similar_sequence(sseqs, query_id, qlen):
    subject_id = sseqs['Hit_id']
    if "BL_ORD_ID" in subject_id:
        subject_id = sseqs['Hit_def'].split(' ')[0]
    if isinstance(sseqs['Hit_hsps']['Hsp'], list):
        coverage = sum([(int(hsp['Hsp_query-to']) -
                         int(hsp['Hsp_query-from']))
                        for hsp in sseqs['Hit_hsps']['Hsp']])
        for hsp in sseqs['Hit_hsps']['Hsp']:
            print_hsp(hsp, query_id, subject_id, coverage, qlen)
    else:
        hsp = sseqs['Hit_hsps']['Hsp']
        coverage = (int(hsp['Hsp_query-to']) -
                    int(hsp['Hsp_query-from']))
        print_hsp(hsp, query_id, subject_id, coverage, qlen)


def test_readxml():
    print()
    inf = "../spec/sample_reports/with_hits_sample/single_iteration.xml"
    tabulardatafromxml(inf)
    inf = "../spec/sample_reports/with_hits_sample/sequenceserver-xml_report.xml"
    tabulardatafromxml(inf)


if __name__ == '__main__':
    import argh

    argh.dispatch_commands([
        tabulardatafromxml
    ])
