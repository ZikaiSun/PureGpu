package additionaluserinterface;

import java.awt.Dimension;
import java.awt.Point;
import java.text.DecimalFormat;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;

/**
 * This class extends JFrame and draw a simple table. All values are in 2D
 * double arrays
 * 
 * @author Daniel Sage, Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
 * 
 */

public class NumericTable extends JFrame {

	private JTable				table;
	private DefaultTableModel	model;

	public NumericTable(String title, String[] headings, Dimension dim) {
		super(title);
		setMinimumSize(dim);
		setSize(dim);
		setPreferredSize(dim);

		JScrollPane pane = new JScrollPane(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS, ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		model = new DefaultTableModel();
		table = new JTable(model);
		for (int i = 0; i < headings.length; i++) {
			model.addColumn(headings[i]);
		}

		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		pane.getViewport().add(table, null);
		add(pane);
	}

	public void setData(double data[][]) {
		int nrow = data.length;
		int ncol = data[0].length;
		String s[] = new String[ncol];
		for (int r = 0; r < nrow; r++) {
			for (int c = 0; c < ncol; c++)
				s[c] = "" + data[r][c];
			model.addRow(s);
		}
	}

	public void setData(double data[][], String[] formats) {
		int nrow = data.length;
		int ncol = data[0].length;
		String s[] = new String[ncol];
		for (int r = 0; r < nrow; r++) {
			for (int c = 0; c < ncol; c++)
				s[c] = (new DecimalFormat(formats[c])).format(data[r][c]);
			model.addRow(s);
		}
	}

	public void setColumnSize(int width[]) {
		for (int i = 0; i < width.length; i++) {
			TableColumn column = table.getColumnModel().getColumn(i);
			column.setPreferredWidth(width[i]);
		}
	}

	public void show(int posx, int posy) {
		pack();
		setLocation(new Point(posx, posy));
		setVisible(true);
	}

}
