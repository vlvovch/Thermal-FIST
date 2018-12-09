#include "tablemodel.h"

#include <QDebug>
#include <QFont>
#include <QDir>
#include <fstream>

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;

const int TableModel::columnNumber = 18;

TableModel::TableModel(QObject *parent)
    :QAbstractTableModel(parent)
{
    RowToParticle.resize(0);
    RowToParticleAll.resize(0);
    RowToParticleStable.resize(0);
    showOnlyStable = false;
}

int TableModel::rowCount(const QModelIndex & /*parent*/) const
{
   return RowToParticle.size();
}

int TableModel::columnCount(const QModelIndex & /*parent*/) const
{
    return columnNumber;
}

QVariant TableModel::data(const QModelIndex &index, int role) const
{

    int row = index.row();
    int col = index.column();
    switch(role){
        case Qt::DisplayRole:
            if (col==0) return QString(model->TPS()->Particles()[RowToParticle[row]].Name().c_str());
            if (col==1) return model->TPS()->Particles()[RowToParticle[row]].PdgId();
            if (col==2) return model->TPS()->Particles()[RowToParticle[row]].Mass();
            if ((col==3 && model->TPS()->Particles()[RowToParticle[row]].IsStable()) ||
                (col==4 && model->TPS()->Particles()[RowToParticle[row]].IsNeutral()!=0)) return "*";
            if (col==5) {
                if (model->TPS()->Particles()[RowToParticle[row]].BaryonCharge()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].BaryonCharge());
                else if (model->TPS()->Particles()[RowToParticle[row]].BaryonCharge()<0) return QString::number(model->TPS()->Particles()[RowToParticle[row]].BaryonCharge());
            }
            if (col==6 && model->TPS()->Particles()[RowToParticle[row]].ElectricCharge()!=0) {
                if (model->TPS()->Particles()[RowToParticle[row]].ElectricCharge()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].ElectricCharge());
                else return model->TPS()->Particles()[RowToParticle[row]].ElectricCharge();
            }
            if (col==7) {
                if (model->TPS()->Particles()[RowToParticle[row]].Strangeness()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].Strangeness());
                else if (model->TPS()->Particles()[RowToParticle[row]].Strangeness()<0) return QString::number(model->TPS()->Particles()[RowToParticle[row]].Strangeness());
                else if (model->TPS()->Particles()[RowToParticle[row]].AbsoluteStrangeness()>1.e-6) return QString::number(model->TPS()->Particles()[RowToParticle[row]].AbsoluteStrangeness()) + " " + tr("(hidden)");
                else if (model->TPS()->Particles()[RowToParticle[row]].Name().substr(0,2)=="K0") return tr("(hidden)");
            }
            if (col==8) {
                if (model->TPS()->Particles()[RowToParticle[row]].Charm()>0) return "+" + QString::number(model->TPS()->Particles()[RowToParticle[row]].Charm());
                else if (model->TPS()->Particles()[RowToParticle[row]].Charm()<0) return QString::number(model->TPS()->Particles()[RowToParticle[row]].Charm());
                else if (model->TPS()->Particles()[RowToParticle[row]].AbsoluteCharm()>1.e-6) return QString::number(model->TPS()->Particles()[RowToParticle[row]].AbsoluteCharm()) + " " + tr("(hidden)");
            }

            if (model->IsCalculated()) {
                if (col==9) return model->GetParticlePrimordialDensity(RowToParticle[row]);
                if (col==10) return model->GetParticlePrimordialDensity(RowToParticle[row]) * model->Volume();
                if (col==11) return model->GetParticleTotalDensity(RowToParticle[row]) * model->Volume();
								if (model->IsFluctuationsCalculated()) {
									if (col == 12) return model->ScaledVariancePrimordial(RowToParticle[row]);
									if (col == 13) return model->ScaledVarianceTotal(RowToParticle[row]);
									if (model->Ensemble() == ThermalModelBase::GCE && model->InteractionModel() == ThermalModelBase::Ideal) {
										if (col == 14) return model->SkewnessPrimordial(RowToParticle[row]);
										if (col == 15) return model->SkewnessTotal(RowToParticle[row]);
										if (col == 16) return model->KurtosisPrimordial(RowToParticle[row]);
										if (col == 17) return model->KurtosisTotal(RowToParticle[row]);
									}
								}
            }
            break;
        case Qt::FontRole:
            if (model->TPS()->Particles()[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
            {
                QFont boldFont;
                boldFont.setBold(true);
                return boldFont;
            }
            break;
        /*case Qt::BackgroundRole:

            if (row == 1 && col == 2)  //change background only for cell(1,2)
            {
                QBrush redBackground(Qt::red);
                return redBackground;
            }
            break;*/
        case Qt::TextAlignmentRole:

            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
            {
                return Qt::AlignCenter;
            }
            break;
        case Qt::CheckStateRole:

            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
            {
                //return Qt::Unchecked;
            }
    }
    return QVariant();
}

QVariant TableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal) {
            switch (section)
            {
                case 1:
                    return tr("PDG ID");
                case 0:
                    return tr("Name");
                case 2:
                    return tr("Mass");
                case 3:
                    return tr("Stable?");
                case 5:
                    return tr("B");
                case 4:
                    return tr("Neutral?");
                case 6:
                    return tr("Q");
                case 7:
                    return tr("S");
                case 8:
                    return tr("C");
                case 9:
                    return tr("Prim. density");
                case 10:
                    return tr("Prim. multiplicity");
                case 11:
                    return tr("Total multiplicity");
                case 12:
                    return tr("Scaled variance prim.");
                case 13:
                    return tr("Scaled variance total");
                case 14:
                    return tr("Skewness prim.");
                case 15:
                    return tr("Skewness total");
                case 16:
                    return tr("Kurtosis prim.");
                case 17:
                    return tr("Kurtosis total");
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}


void TableModel::reset() {
    beginResetModel();
    RowToParticle.resize(0);
    RowToParticleAll.resize(0);
    RowToParticleStable.resize(0);
    for(int i=0;i<model->TPS()->Particles().size();++i) {
        RowToParticleAll.push_back(i);
        if (model->TPS()->Particles()[i].IsStable()) RowToParticleStable.push_back(i);
    }
    if (showOnlyStable) RowToParticle = RowToParticleStable;
    else RowToParticle = RowToParticleAll;
    endResetModel();
}

void TableModel::setModel(ThermalModelBase *mod) {
    beginResetModel();
    model = mod;
    RowToParticle.resize(0);
    RowToParticleAll.resize(0);
    RowToParticleStable.resize(0);
    for(int i=0;i<model->TPS()->Particles().size();++i) {
        RowToParticleAll.push_back(i);
        if (model->TPS()->Particles()[i].IsStable()) RowToParticleStable.push_back(i);
    }
    if (showOnlyStable) RowToParticle = RowToParticleStable;
    else RowToParticle = RowToParticleAll;
    endResetModel();
}


DecayTableModel::DecayTableModel(QObject *parent, ThermalParticle *Particle, ThermalParticleSystem *TPS)
    :QAbstractTableModel(parent), fParticle(Particle), fTPS(TPS)
{
    columnNumber = 3;
    bratioSum = 0.;

    for(unsigned int i = 0; i < Particle->Decays().size(); i++) {
        if (Particle->Decays()[i].mDaughters.size()+2>columnNumber) columnNumber = 2 + Particle->Decays()[i].mDaughters.size();
        bratioSum += Particle->Decays()[i].mBratio;
    }
}


int DecayTableModel::rowCount(const QModelIndex & /*parent*/) const
{
   return fParticle->Decays().size();
}

int DecayTableModel::columnCount(const QModelIndex & /*parent*/) const
{
    return columnNumber;
}

QVariant DecayTableModel::data(const QModelIndex &index, int role) const
{
    int row = index.row();
    int col = index.column();
    ParticleDecay decay;
    switch(role){
        case Qt::DisplayRole:
            if (row>=fParticle->Decays().size()) return QVariant();
            decay = fParticle->Decays()[row];
            if (col==0) return decay.mBratio * 100.;
            else if (col==1) return decay.mBratio / bratioSum * 100.;
            else if ((col-2)<decay.mDaughters.size()){
                int PDGID = decay.mDaughters[col-2];
                if (fTPS->PdgToId(PDGID)==-1) return QString("%1   %2").arg(QString::number(PDGID), "???");
                else return QString("%1   %2").arg(QString::number(PDGID), QString(fTPS->ParticleByPDG(PDGID).Name().c_str()));
            }
            else return QVariant();
            break;
        case Qt::EditRole:
            if (row>=fParticle->Decays().size()) return QVariant();
            decay = fParticle->Decays()[row];
            if (col==0) return decay.mBratio * 100.;
            else if (col==1) return decay.mBratio / bratioSum * 100.;
            else if ((col-2)<decay.mDaughters.size()){
                return decay.mDaughters[col-2];
            }
            else return 0;
            break;
        /*case Qt::FontRole:
            if (TPS.fParticles[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
            {
                QFont boldFont;
                boldFont.setBold(true);
                return boldFont;
            }
            break;
        case Qt::BackgroundRole:

            if (row == 1 && col == 2)  //change background only for cell(1,2)
            {
                QBrush redBackground(Qt::red);
                return redBackground;
            }
            break;
        case Qt::TextAlignmentRole:

            if (col>=3 && col<=6) //add a checkbox to cell(1,0)
            {
                return Qt::AlignCenter;
            }
            break;
        case Qt::CheckStateRole:

            if (col>=3 && col<=6) //add a checkbox to cell(1,0)
            {
                //return Qt::Unchecked;
            }*/
    }
    return QVariant();
}

QVariant DecayTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal) {
            //switch (section)
            {
                if (section==0) return tr("Branching ratio (%)");
                if (section==1) return tr("Normalized (%)");
                if (section>=2) return QString(tr("Daughter %1")).arg(section-1);
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}

bool DecayTableModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
    int row = index.row();
    int col = index.column();
    //qDebug() << value.toString();
    if (role==Qt::EditRole) {
        if (col==0) {
            fParticle->Decays()[row].mBratio = value.toDouble() / 100.;
            return true;
        }
        else if (col!=1 && (col-2)<fParticle->Decays()[row].mDaughters.size()) {
            //qDebug() << col - 2 ;// << " " << distance(fParticle->Decays()[row].mDaughters.begin(), fParticle->Decays()[row].mDaughters.begin() + col - 2);
            if (value.toInt()!=0)
                fParticle->Decays()[row].mDaughters[col-2] = value.toInt();
            else fParticle->Decays()[row].mDaughters.erase(fParticle->Decays()[row].mDaughters.begin() + (int)(col - 2));
            return true;
        }
        else if (col!=1 && value.toInt()!=0) {
            fParticle->Decays()[row].mDaughters.push_back(value.toInt());
            return true;
        }
    }
    else return false;
}


QuantitiesModel::QuantitiesModel(QObject *parent, std::vector<FittedQuantity>  *quants, ThermalModelBase *modelop)
    :quantities(quants), model(modelop)
{
}


int QuantitiesModel::rowCount(const QModelIndex & /*parent*/) const
{
   return quantities->size();
}

int QuantitiesModel::columnCount(const QModelIndex & /*parent*/) const
{
    return 8;
}

QVariant QuantitiesModel::data(const QModelIndex &index, int role) const
{
    int row = index.row();
    int col = index.column();
    FittedQuantity quantity;
    switch(role){
        case Qt::DisplayRole:
            if (row>=quantities->size()) return QVariant();
            quantity = (*quantities)[row];
            if (col==0) {
                if (quantity.type == FittedQuantity::Multiplicity) {
                    return QString::fromStdString(model->TPS()->GetNameFromPDG(quantity.mult.fPDGID));
                }
                else {
                    QString name1, name2;
                    name1 = QString::fromStdString(model->TPS()->GetNameFromPDG(quantity.ratio.PDGID1));
                    name2 = QString::fromStdString(model->TPS()->GetNameFromPDG(quantity.ratio.PDGID2));

                    return name1 + "/" + name2;
                }
            }
            else if (col==2) {
                if (quantity.type == FittedQuantity::Multiplicity) return quantity.mult.fValue;
                else return quantity.ratio.fValue;
            }
            else if (col==3) {
                if (quantity.type == FittedQuantity::Multiplicity) return quantity.mult.fError;
                else return quantity.ratio.fError;
            }
            else if (col==4 && model->IsCalculated()) {
                if (quantity.type == FittedQuantity::Multiplicity) {
                    double val = model->GetDensity(quantity.mult.fPDGID, quantity.mult.fFeedDown);

                    return val * model->Parameters().V;
                }
                else {
                    double val1 = model->GetDensity(quantity.ratio.PDGID1, quantity.ratio.fFeedDown1);
                    double val2 = model->GetDensity(quantity.ratio.PDGID2, quantity.ratio.fFeedDown2);

                    return val1 / val2;
                }
            }
            else if (col==5 && model->IsCalculated()) {
                if (quantity.type == FittedQuantity::Multiplicity) {

                    double val = model->GetDensity(quantity.mult.fPDGID, quantity.mult.fFeedDown);
                    val = val * model->Parameters().V;

                    double valexp = quantity.mult.fValue;
                    double err = quantity.mult.fError;

                    return (val - valexp) / err;
                }
                else {
                    double val1 = model->GetDensity(quantity.ratio.PDGID1, quantity.ratio.fFeedDown1);
                    double val2 = model->GetDensity(quantity.ratio.PDGID2, quantity.ratio.fFeedDown2);

                    double valexp = quantity.ratio.fValue;
                    double err = quantity.ratio.fError;

                    return (val1 / val2 - valexp) / err;
                }
            }
						else if (col == 6 && model->IsCalculated()) {
							if (quantity.type == FittedQuantity::Multiplicity) {
								double val = model->GetDensity(quantity.mult.fPDGID, quantity.mult.fFeedDown);
								val = val * model->Parameters().V;

								double valexp = quantity.mult.fValue;
								double err = quantity.mult.fError;

								QString str1 = QString::number(valexp / val, 'f', 3);
								QString str2 = QString::number(err / val, 'f', 3);

								return str1 + " ± " + str2;
							}
							else {
								double val1 = model->GetDensity(quantity.ratio.PDGID1, quantity.ratio.fFeedDown1);
								double val2 = model->GetDensity(quantity.ratio.PDGID2, quantity.ratio.fFeedDown2);

								double valexp = quantity.ratio.fValue;
								double err = quantity.ratio.fError;

								QString str1 = QString::number(valexp / (val1/val2), 'f', 3);
								QString str2 = QString::number(err / (val1 / val2), 'f', 3);

								return str1 + " ± " + str2;
							}
						}
            else if (col==7) {
                if (quantity.type == FittedQuantity::Multiplicity) {
                    if (quantity.mult.fFeedDown==0) return QString(tr("Primordial"));
                    if (quantity.mult.fFeedDown==1) return QString(tr("Strong decays"));
                    if (quantity.mult.fFeedDown==2) return QString(tr("Strong+weak decays"));
                }
                else {
                    QString str1, str2;
                    if (quantity.ratio.fFeedDown1==0) str1 = QString(tr("Primordial"));
                    if (quantity.ratio.fFeedDown1==1) str1 = QString(tr("Strong decays"));
                    if (quantity.ratio.fFeedDown1==2) str1 = QString(tr("Strong+weak decays"));
                    if (quantity.ratio.fFeedDown2==0) str2 = QString(tr("Primordial"));
                    if (quantity.ratio.fFeedDown2==1) str2 = QString(tr("Strong decays"));
                    if (quantity.ratio.fFeedDown2==2) str2 = QString(tr("Strong+weak decays"));
                    return str1 + "/" + str2;
                }
                return QVariant();
            }
            else return QVariant();
            break;
        /*case Qt::EditRole:
            if (row>=fParticle->Decays().size()) return QVariant();
            decay = fParticle->Decays()[row];
            if (col==0) return decay.mBratio * 100.;
            else if (col==1) return decay.mBratio / bratioSum * 100.;
            else if ((col-2)<decay.mDaughters.size()){
                return decay.mDaughters[col-2];
            }
            else return 0;
            break;
        case Qt::FontRole:
            if (TPS.fParticles[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
            {
                QFont boldFont;
                boldFont.setBold(true);
                return boldFont;
            }
            break;
        case Qt::BackgroundRole:

            if (row == 1 && col == 2)  //change background only for cell(1,2)
            {
                QBrush redBackground(Qt::red);
                return redBackground;
            }
            break;*/
        case Qt::TextAlignmentRole:

            if (col == 1) //add a checkbox to cell(1,0)
            {
							return Qt::AlignCenter;
            }
            break;
        case Qt::CheckStateRole:
					if (row >= quantities->size()) return QVariant();
					quantity = (*quantities)[row];
          if (col == 1) //add a checkbox to cell(1,0)
          {
              if (quantity.toFit)
								return Qt::Checked;
							else
								return Qt::Unchecked;
          }
    }
    return QVariant();
}

QVariant QuantitiesModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal) {
            //switch (section)
            {
                if (section == 0) return tr("Name");
								if (section == 1) return tr("Fit?");
                if (section == 2) return tr("Exp. value");
                if (section == 3) return tr("Exp. error");
                if (section == 4) return tr("Model value");
                if (section == 5) return tr("Deviation");
                if (section == 6) return tr("Data/Model");
								if (section == 7) return tr("Feeddown");
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}

Qt::ItemFlags QuantitiesModel::flags(const QModelIndex &index) const
{
	Qt::ItemFlags flags = QAbstractItemModel::flags(index);
	if (index.column() == 1)
	{
		flags |= Qt::ItemIsUserCheckable;
		flags |= Qt::ItemIsSelectable;
	}
	return flags;
}

bool QuantitiesModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
	int row = index.row();
	int col = index.column();
	if (row >= quantities->size())
		return false;
	FittedQuantity& quantity = (*quantities)[row];
	if (role == Qt::CheckStateRole) {
		quantity.toFit = value.toBool();
	}
	emit dataChanged(index, index);
	return true;
}

void QuantitiesModel::setModel(ThermalModelBase *mod) {
    beginResetModel();
    model = mod;
    endResetModel();
}

void QuantitiesModel::setQuantities(std::vector<FittedQuantity>  *quants) {
    beginResetModel();
    quantities = quants;
    endResetModel();
}



const int SpectraModel::columnNumber = 9;

SpectraModel::SpectraModel(QObject *parent)
    :QAbstractTableModel(parent)
{
    RowToParticle.resize(0);
}

int SpectraModel::rowCount(const QModelIndex & /*parent*/) const
{
    if (RowToParticle.size()>0) return RowToParticle.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + spectra->fNegativeCharges.size();
    return RowToParticle.size();
}

int SpectraModel::columnCount(const QModelIndex & /*parent*/) const
{
    return columnNumber;
}

QVariant SpectraModel::data(const QModelIndex &index, int role) const
{
    int row = index.row();
    int col = index.column();
    switch(role){
        case Qt::DisplayRole:
            if (col==0) {
                if (row<spectra->fNames.size()) return QString(spectra->fNames[RowToParticle[row]].c_str());
                if (row<spectra->fNames.size() + spectra->fNetParticles.size())
                    return spectra->fNetParticles[row-spectra->fNames.size()].GetName().c_str();
                if (row==spectra->fNames.size() + spectra->fNetParticles.size()) return QString("net-baryon");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + 1) return QString("net-charge");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + 2) return QString("net-strangeness");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size()) return QString("baryon hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + 1) return QString("charge hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + 2) return QString("strange hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size()) return QString("baryon+ hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + 1) return QString("charge+ hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + 2) return QString("strange+ hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size()) return QString("baryon- hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + 1) return QString("charge- hadrons");
                if (row==spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + 2) return QString("strange- hadrons");
            }
            if (col==1) {
                if (row<spectra->fParticles.size()) return spectra->fParticles[RowToParticle[row]].GetPDGID();
                else return "--";
            }
            if (col==2) {
                if (row<spectra->fMasses.size()) return spectra->fMasses[RowToParticle[row]];
                else return "--";
            }
            {
                if (row<spectra->fParticles.size()) {
                    if (spectra->fParticles[RowToParticle[row]].GetEvents()==0) return QVariant();
                    if (col==3) return QString::number(spectra->fParticles[RowToParticle[row]].GetMean(),'f',3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetMeanError(),'f',3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
                    if (col==4) return QString::number(spectra->fParticles[RowToParticle[row]].GetVariance(),'f',3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetVarianceError(),'f',3);
                    if (col==5) return QString::number(spectra->fParticles[RowToParticle[row]].GetScaledVariance(),'f',3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetScaledVarianceError(),'f',3);
                    if (col==6) return QString::number(spectra->fParticles[RowToParticle[row]].GetSkewness(),'f',3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetSkewnessError(),'f',3);
                    if (col==7) return QString::number(spectra->fParticles[RowToParticle[row]].GetKurtosis(),'f',3) + " ± " + QString::number(spectra->fParticles[RowToParticle[row]].GetKurtosisError(),'f',3);
                }
                else if (row<spectra->fNames.size() + spectra->fNetParticles.size()) {
                    int tind = row - spectra->fParticles.size();
                    if (spectra->fNetParticles[tind].n==0.) return QVariant();
                    if (col==3) return QString::number(spectra->fNetParticles[tind].GetMean(),'f',3) + " ± " + QString::number(spectra->fNetParticles[tind].GetMeanError(),'f',3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
                    if (col==4) return QString::number(spectra->fNetParticles[tind].GetVariance(),'f',3) + " ± " + QString::number(spectra->fNetParticles[tind].GetVarianceError(),'f',3);
                    if (col==5) return QString::number(spectra->fNetParticles[tind].GetScaledVariance(),'f',3) + " ± " + QString::number(spectra->fNetParticles[tind].GetScaledVarianceError(),'f',3);
                    if (col==6) return QString::number(spectra->fNetParticles[tind].GetSkewness(),'f',3) + " ± " + QString::number(spectra->fNetParticles[tind].GetSkewnessError(),'f',3);
                    if (col==7) return QString::number(spectra->fNetParticles[tind].GetKurtosis(),'f',3) + " ± " + QString::number(spectra->fNetParticles[tind].GetKurtosisError(),'f',3);
                }
                else if (row<spectra->fNames.size() + spectra->fNetParticles.size() + 3) {
                    int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size();
                    if (spectra->fNetCharges[tind].n==0.) return QVariant();
                    if (col==3) return QString::number(spectra->fNetCharges[tind].GetMean(),'f',3) + " ± " + QString::number(spectra->fNetCharges[tind].GetMeanError(),'f',3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
                    if (col==4) return QString::number(spectra->fNetCharges[tind].GetVariance(),'f',3) + " ± " + QString::number(spectra->fNetCharges[tind].GetVarianceError(),'f',3);
                    if (col==5) return QString::number(spectra->fNetCharges[tind].GetScaledVariance(),'f',3) + " ± " + QString::number(spectra->fNetCharges[tind].GetScaledVarianceError(),'f',3);
                    if (col==6) return QString::number(spectra->fNetCharges[tind].GetSkewness(),'f',3) + " ± " + QString::number(spectra->fNetCharges[tind].GetSkewnessError(),'f',3);
                    if (col==7) return QString::number(spectra->fNetCharges[tind].GetKurtosis(),'f',3) + " ± " + QString::number(spectra->fNetCharges[tind].GetKurtosisError(),'f',3);
                }
                else if (row<spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + 3) {
                    int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size() - spectra->fNetCharges.size();
                    if (spectra->fTotalCharges[tind].n==0.) return QVariant();
                    if (col==3) return QString::number(spectra->fTotalCharges[tind].GetMean(),'f',3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetMeanError(),'f',3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
                    if (col==4) return QString::number(spectra->fTotalCharges[tind].GetVariance(),'f',3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetVarianceError(),'f',3);
                    if (col==5) return QString::number(spectra->fTotalCharges[tind].GetScaledVariance(),'f',3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetScaledVarianceError(),'f',3);
                    if (col==6) return QString::number(spectra->fTotalCharges[tind].GetSkewness(),'f',3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetSkewnessError(),'f',3);
                    if (col==7) return QString::number(spectra->fTotalCharges[tind].GetKurtosis(),'f',3) + " ± " + QString::number(spectra->fTotalCharges[tind].GetKurtosisError(),'f',3);
                }
                else if (row<spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + 3) {
                    int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size() - spectra->fNetCharges.size() - spectra->fTotalCharges.size();
                    if (spectra->fPositiveCharges[tind].n==0.) return QVariant();
                    if (col==3) return QString::number(spectra->fPositiveCharges[tind].GetMean(),'f',3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetMeanError(),'f',3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
                    if (col==4) return QString::number(spectra->fPositiveCharges[tind].GetVariance(),'f',3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetVarianceError(),'f',3);
                    if (col==5) return QString::number(spectra->fPositiveCharges[tind].GetScaledVariance(),'f',3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetScaledVarianceError(),'f',3);
                    if (col==6) return QString::number(spectra->fPositiveCharges[tind].GetSkewness(),'f',3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetSkewnessError(),'f',3);
                    if (col==7) return QString::number(spectra->fPositiveCharges[tind].GetKurtosis(),'f',3) + " ± " + QString::number(spectra->fPositiveCharges[tind].GetKurtosisError(),'f',3);
                }
                else if (row<spectra->fNames.size() + spectra->fNetParticles.size() + spectra->fNetCharges.size() + spectra->fTotalCharges.size() + spectra->fPositiveCharges.size() + 3) {
                    int tind = row - spectra->fParticles.size() - spectra->fNetParticles.size() - spectra->fNetCharges.size() - spectra->fTotalCharges.size() - spectra->fPositiveCharges.size();
                    if (spectra->fNegativeCharges[tind].n==0.) return QVariant();
                    if (col==3) return QString::number(spectra->fNegativeCharges[tind].GetMean(),'f',3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetMeanError(),'f',3);//QString("%1 ± %2").arg(spectra->fParticles[RowToParticle[row]].GetMean()).arg(spectra->fParticles[RowToParticle[row]].GetMeanError());
                    if (col==4) return QString::number(spectra->fNegativeCharges[tind].GetVariance(),'f',3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetVarianceError(),'f',3);
                    if (col==5) return QString::number(spectra->fNegativeCharges[tind].GetScaledVariance(),'f',3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetScaledVarianceError(),'f',3);
                    if (col==6) return QString::number(spectra->fNegativeCharges[tind].GetSkewness(),'f',3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetSkewnessError(),'f',3);
                    if (col==7) return QString::number(spectra->fNegativeCharges[tind].GetKurtosis(),'f',3) + " ± " + QString::number(spectra->fNegativeCharges[tind].GetKurtosisError(),'f',3);
                }
            }
            if (col==8 && row<spectra->fParticles.size() && spectra->fParticles[RowToParticle[row]].GetAcceptance()) return "*";
            break;
        case Qt::FontRole:
//            if (model->TPS()->Particles()[RowToParticle[row]].IsStable()) //change font only for cell(0,0)
//            {
//                QFont boldFont;
//                boldFont.setBold(true);
//                return boldFont;
//            }
            break;
        /*case Qt::BackgroundRole:

            if (row == 1 && col == 2)  //change background only for cell(1,2)
            {
                QBrush redBackground(Qt::red);
                return redBackground;
            }
            break;*/
        case Qt::TextAlignmentRole:

//            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
//            {
//                return Qt::AlignCenter;
//            }
              if (col==8) return Qt::AlignCenter;
            break;
        case Qt::CheckStateRole:

//            if (col>=3 && col<=8) //add a checkbox to cell(1,0)
//            {
//                //return Qt::Unchecked;
//            }
        break;
    }
    return QVariant();
}

QVariant SpectraModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal) {
            switch (section)
            {
                case 1:
                    return tr("PDG ID");
                case 0:
                    return tr("Name");
                case 2:
                    return tr("Mass");
                case 3:
                    return tr("Multiplicity");
                case 4:
                    return tr("Variance");
                case 5:
                    return tr("Scaled variance");
                case 6:
                    return tr("Skewness");
                case 7:
                    return tr("Kurtosis");
                case 8:
                    return tr("Acceptance");
            }
        }
        else if (orientation == Qt::Vertical) return section + 1;
    }
    return QVariant();
}


void SpectraModel::reset() {
    beginResetModel();
    //TPS = TPScopy;
    RowToParticle.resize(0);
    for(int i=0;i<spectra->fParticles.size();++i) {
        RowToParticle.push_back(i);
    }
//    if (showOnlyStable) RowToParticle = RowToParticleStable;
//    else RowToParticle = RowToParticleAll;
    endResetModel();
    //qDebug() << RowToParticle.size();
}

void SpectraModel::setSpectra(ParticlesSpectra *spec) {
    beginResetModel();
    spectra = spec;
    RowToParticle.resize(0);
    if (spectra!=NULL) {
        for(int i=0;i<spectra->fParticles.size();++i) {
            RowToParticle.push_back(i);
        }
    }
    endResetModel();
}
